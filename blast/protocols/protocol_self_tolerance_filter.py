# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Enzo Sierra (enzogael57@gmail.com)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This protocol re-checks self-tolerance (human-proteome homology) at the
FINAL candidate's own length, rather than at the length of whatever
upstream parent region it was selected from.

Rationale: an earlier BLASTp self-tolerance check run against a longer
PARENT region (e.g. a fused B-cell union region, up to ~20-70aa) can miss
a short (8-15aa) dangerous human-homologous motif buried inside it,
because the minimum-query-coverage filter is computed against the
parent's full length, not the short motif's own length. This protocol
re-runs the same BLASTp self-tolerance check on the actual FINAL
candidate sequence -- coverage is now measured against ITS OWN length, at
the scale where it actually matters.

Reusable across any candidate class (B-cell/HTL/CTL): unlike the crossref
protocols (ProtLANLCATNAPCrossref/ProtIEDBCrossref), which only annotate,
this one FILTERS (drops candidates whose max human-proteome identity
exceeds the threshold). Meant to be inserted at more than one point of a
workflow -- once on the longer parent regions, again on the final
per-class candidate sequences -- to close the coverage gap described
above.
"""

import os
from pathlib import Path
from typing import List

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, String
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from blast import Plugin as blastPlugin
from blast.constants import BLAST_HUMAN_DB_PATH
from blast.utils.blastp_batch import maxIdentityByQuery, runBatchedBlastp

DEFAULT_IDENTITY_THRESHOLD = 75.0


class ProtSelfToleranceFilter(EMProtocol):
    """
    AI Generated:

    Re-checks every input ROI's peptide against a local human-proteome
    BLAST database, at the peptide's OWN length (not an upstream parent
    region's length -- see module docstring), and keeps only ROIs whose
    maximum human-proteome identity does NOT exceed 'identityThreshold'.

    Output
    ------
    outputROIs: the SURVIVING subset of the input SetOfSequenceROIs (this
    protocol filters, it does not just annotate), each kept ROI annotated
    with '_maxHumanPident' (float, 0.0 if no qualifying hit) and
    '_blastVerdict' (str, always 'Segura' for a surviving ROI -- kept for
    traceability, consistent with the verdict column of other BLAST-based
    protocols in this plugin).
    Full per-candidate detail (survivors and rejects alike) persisted to
    'extra/self_tolerance_check.csv'.
    """

    _label = 'self tolerance filter'

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Final candidate peptides to re-check at their own length (e.g. '
                            'HTL/CTL cores, or a B-cell candidate already trimmed/padded).')
        form.addParam('humanDbPath', params.StringParam, default='',
                       label='Human proteome BLAST DB prefix (blank = scipion.conf default): ',
                       help=f'Prefix (no extension) of a local, already-indexed BLASTp protein '
                            f'database of the human proteome. Blank uses the {BLAST_HUMAN_DB_PATH} '
                            f'variable from scipion.conf.')
        form.addParam('identityThreshold', params.FloatParam, default=DEFAULT_IDENTITY_THRESHOLD,
                       label='Max. human identity %% (exclusive): ',
                       help='A candidate is discarded if its best qualifying BLASTp hit against '
                            'the human proteome exceeds this %% identity.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.filterStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getResultPath(self):
        return self._getExtraPath('self_tolerance_check.csv')

    def _getRois(self) -> List:
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def _getDbPath(self) -> str:
        return self.humanDbPath.get().strip() or blastPlugin.getVar(BLAST_HUMAN_DB_PATH)

    def filterStep(self):
        rois = self._getRois()
        sequences = pd.Series([roi.getROISequence() for roi in rois])
        if sequences.empty:
            pd.DataFrame(columns=['sequence', 'max_pident', 'status']).to_csv(self._getResultPath(), index=False)
            return

        hits = runBatchedBlastp(self, blastPlugin, sequences, self._getDbPath())
        maxPident = maxIdentityByQuery(hits, sequences.str.len(), minQueryCoverage=0.9)

        threshold = self.identityThreshold.get()
        result = pd.DataFrame({'sequence': sequences})
        result['max_pident'] = [float(maxPident.get(f'peptide_{idx}', 0.0)) for idx in sequences.index]
        result['status'] = result['max_pident'].apply(lambda p: 'Autoinmunidad' if p > threshold else 'Segura')
        result.to_csv(self._getResultPath(), index=False)

    def createOutputStep(self):
        rois = self._getRois()
        resultDf = pd.read_csv(self._getResultPath()) if os.path.isfile(self._getResultPath()) else pd.DataFrame()

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for i, roi in enumerate(rois):
            if resultDf.empty or i >= len(resultDf):
                continue
            row = resultDf.iloc[i]
            if row['status'] != 'Segura':
                continue
            roi._maxHumanPident = Float(float(row['max_pident']))
            roi._blastVerdict = String('Segura')
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        errors = []
        dbPath = self._getDbPath()
        if not dbPath:
            errors.append(
                f'No human proteome BLAST DB configured: set humanDbPath or the '
                f'{BLAST_HUMAN_DB_PATH} variable in scipion.conf.'
            )
        elif not Path(f'{dbPath}.phr').is_file():
            errors.append(f"Human proteome BLAST DB not found at '{dbPath}' (expected '{dbPath}.phr').")
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            resultDf = pd.read_csv(self._getResultPath()) if os.path.isfile(self._getResultPath()) else pd.DataFrame()
            outROIs = getattr(self, 'outputROIs', None)
            nSurvivors = len(outROIs) if outROIs is not None else 0
            nTotal = len(resultDf)
            summary.append(f'{nSurvivors}/{nTotal} candidate(s) survive the self-tolerance re-check '
                            f'at their own length.')
        return summary

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
This protocol annotates each input ROI with how conserved it is across a
user-provided reference panel of strain/clade/variant sequences of the
SAME pathogen under analysis.

Ported from the standalone B-Cell-Epitope-Prediction repo's
src/engines/conservation_engine.py (Fase 6b): a poorly-conserved epitope
protects only against the exact strain it was designed against;
prioritizing conserved epitopes broadens
coverage and puts more evolutionary pressure on the pathogen (conserved
regions tend to be functionally constrained -- mutating them has a real
fitness cost). Unlike ProtSelfToleranceFilter's fixed human proteome
database, there is no universal conservation panel: it is specific to
the strains/variants of whichever pathogen is under analysis in a given
run, provided fresh each time (a SetOfSequences the user assembled/
imported), not a fixed scipion.conf variable.

Metric is BREADTH, not best-hit: counts how many DISTINCT panel sequences
a candidate matches above the identity/coverage thresholds, not just the
single closest one. A candidate identical to one rare panel strain is not
"conserved" in the sense that matters here (protecting against MANY
circulating variants); one with an acceptable match against 80%% of the
panel is, even if none of those individual matches is 100%% identity.

Purely informative: never filters (same treatment as
ProtLANLCATNAPCrossref/ProtIEDBCrossref), just annotates.
"""

import hashlib
import os
from pathlib import Path

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwchem.utils.utilsFasta import fastFastaExport
from pwem.protocols import EMProtocol
from pyworkflow.object import Float, Integer
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from blast import Plugin as blastPlugin
from blast.utils.blastp_batch import panel_breadth_by_query, run_batched_blastp

DEFAULT_IDENTITY_THRESHOLD = 90.0
DEFAULT_MIN_QUERY_COVERAGE = 0.9


class ProtBLASTPanelConservation(EMProtocol):
    """
    AI Generated:

    Annotates (does NOT filter) every input ROI with its conservation
    breadth against a reference panel of same-pathogen strain/clade/
    variant sequences: the panel is BLASTp-indexed once (cached by content
    hash across re-runs with the same panel) and each ROI's peptide is
    searched against it in length-tiered batches.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, annotated with
    '_conservationPanelMatches' (int, distinct panel sequences matched),
    '_conservationPanelTotal' (int, panel size), '_conservationPct'
    (float, matches/total * 100). Full detail persisted to
    'extra/conservation.csv'.
    """

    _label = 'blast panel conservation'

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to annotate with panel conservation breadth.')
        form.addParam('panelSequences', params.PointerParam, pointerClass='SetOfSequences',
                       label='Reference panel (strains/clades/variants): ',
                       help='Unindexed reference sequences of the SAME pathogen under analysis '
                            '(other strains/clades/variants) -- NOT a fixed database, specific '
                            'to this run.')
        form.addParam('identityThreshold', params.FloatParam, default=DEFAULT_IDENTITY_THRESHOLD,
                       label='Min. identity %% to count a panel match: ',
                       help='A panel sequence counts as matched if a qualifying BLASTp hit '
                            'reaches at least this %% identity.')
        form.addParam('minQueryCoverage', params.FloatParam, default=DEFAULT_MIN_QUERY_COVERAGE,
                       label='Min. query coverage (0-1): ',
                       help='Minimum fraction of the candidate\'s own length a BLASTp alignment '
                            'must cover to count -- without this, a tiny fragment identical by '
                            'chance would count as a real match.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.indexPanelStep)
        self._insertFunctionStep(self.annotateStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getResultPath(self):
        return self._getExtraPath('conservation.csv')

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def indexPanelStep(self):
        # ABSOLUTE path: runBLAST's runJob is invoked with cwd=cacheDir below, so a relative
        # path here would resolve against the WRONG directory (same class of bug already found
        # and fixed in AlgPred2/ToxinPred2/IApred in the standalone project).
        panelFasta = os.path.abspath(self._getExtraPath('panel.fasta'))
        fastFastaExport(self.panelSequences.get(), panelFasta)

        self._nPanelTotal = sum(1 for _ in self.panelSequences.get())
        contentHash = hashlib.sha256(Path(panelFasta).read_bytes()).hexdigest()[:16]

        cacheDir = Path(blastPlugin.getVar('BLAST_HOME') or '.').parent / 'conservation_panel_cache' / contentHash
        cacheDir.mkdir(parents=True, exist_ok=True)
        self._dbPrefix = str(cacheDir / 'panel_db')

        if not Path(f'{self._dbPrefix}.phr').is_file():
            args = f'-in {panelFasta} -parse_seqids -dbtype prot -out {self._dbPrefix}'
            blastPlugin.runBLAST(self, 'makeblastdb', args, cwd=str(cacheDir))

    def annotateStep(self):
        rois = self._getRois()
        sequences = pd.Series([roi.getROISequence() for roi in rois])
        result = pd.DataFrame({'sequence': sequences})

        if sequences.empty:
            result = result.assign(n_panel_matches=[], n_panel_total=[], conservation_pct=[])
        else:
            hits = run_batched_blastp(self, blastPlugin, sequences, self._dbPrefix)
            breadth = panel_breadth_by_query(
                hits, sequences.str.len(), self.identityThreshold.get(), self.minQueryCoverage.get(),
            )
            result['n_panel_matches'] = [int(breadth.get(f'peptide_{idx}', 0)) for idx in sequences.index]
            result['n_panel_total'] = self._nPanelTotal
            result['conservation_pct'] = (result['n_panel_matches'] / self._nPanelTotal * 100.0).round(2)

        result.to_csv(self._getResultPath(), index=False)

    def createOutputStep(self):
        rois = self._getRois()
        resultDf = pd.read_csv(self._getResultPath()) if os.path.isfile(self._getResultPath()) else pd.DataFrame()

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for i, roi in enumerate(rois):
            if not resultDf.empty and i < len(resultDf):
                row = resultDf.iloc[i]
                roi._conservationPanelMatches = Integer(int(row['n_panel_matches']))
                roi._conservationPanelTotal = Integer(int(row['n_panel_total']))
                roi._conservationPct = Float(float(row['conservation_pct']))
            else:
                roi._conservationPanelMatches = Integer(0)
                roi._conservationPanelTotal = Integer(0)
                roi._conservationPct = Float(0.0)
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)
            self._defineSourceRelation(self.panelSequences, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        errors = []
        if self.panelSequences.get() is not None and sum(1 for _ in self.panelSequences.get()) == 0:
            errors.append('The reference panel is empty.')
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            resultDf = pd.read_csv(self._getResultPath()) if os.path.isfile(self._getResultPath()) else pd.DataFrame()
            if not resultDf.empty:
                summary.append(f'Reference panel: {int(resultDf["n_panel_total"].iloc[0])} sequence(s). '
                                f'Mean conservation: {resultDf["conservation_pct"].mean():.2f}%.')
        return summary

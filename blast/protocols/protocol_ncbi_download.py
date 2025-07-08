# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os, json
from Bio import Entrez

from pwem.protocols import EMProtocol
from pwem.objects import Sequence, SetOfSequences
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pyworkflow.protocol.params import TextParam, StringParam, EnumParam, STEPS_PARALLEL, IntParam, LabelParam
from pyworkflow import BETA

IDS, KEYS = 0, 1

class ProtChemNCBIDownload(EMProtocol):
    """Download the Fasta file(s) from NCBI databases

User IA Manual: NcbiDownload Protocol

The NcbiDownload protocol allows users to retrieve biological sequences directly
from the NCBI databases using accession numbers or search terms. It provides a
simple and reproducible way to access curated nucleotide or protein sequences
from GenBank, RefSeq, or UniProt without leaving the Scipion-Chem environment.

To begin, the user must provide a list of identifiers or a search query that
matches entries in the NCBI database. These identifiers can refer to accession
codes, gene symbols, protein names, or any valid term accepted by the NCBI
Entrez system. The user must also select the type of sequence to be retrieved,
choosing between nucleotide and protein formats. This selection determines which
databases will be queried and how the results will be parsed.

The protocol offers the option to download either a single sequence per query or
to collect all matching records, depending on the specificity of the search. In
addition, the user can limit the number of returned sequences to avoid retrieving
large datasets unintentionally. The output format is typically FASTA, and each
entry includes standard annotations such as the sequence identifier, description,
and source organism.

Downloaded sequences are saved locally and registered as Scipion objects, making
them immediately available for use in downstream protocols. These may include
sequence alignment, BLAST database creation, modeling, or mutation analysis.
A log of the download process is also generated, including information about
successful retrievals and any queries that failed to return valid results.

In summary, the NcbiDownload protocol provides a streamlined method for accessing
sequence data from NCBI within the Scipion-Chem framework. It enables automated
and traceable integration of external biological information into
structure-based and bioinformatics workflows.
    
    """
    _label = 'NCBI download'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input type')
        group.addParam('searchMode', EnumParam, default=0, choices=['ID', 'Keyword'], label='Search mode: ',
                      help='Whether to fecth the NCBI entry for an specific ID or to use esearch to find entries '
                           'for a specific keyword')

        group.addParam('dbType', EnumParam, default=0, choices=['Protein', 'Nucleotide', 'PCCompound'],
                       display=EnumParam.DISPLAY_HLIST, label='Type of database to query on: ',
                       help='Type of data to download from NCBI databases. Compounds are fetched from PubChem.')

        group = form.addGroup('Input listing')
        group.addParam('inputID', StringParam, label='NCBI ID / keyword: ', help="NCBI ID / keyword for the query")
        group.addParam('maxEntries', IntParam, label='Maximum number of entries: ', condition="searchMode==1",
                       default=5, help="maximum number of entries to take into account in the search")
        group.addParam('addEntry', LabelParam, label='Add ID / keyword: ', help='Add ID / keyword to the list')

        group.addParam('listIDs', TextParam, width=60, label='List of IDs / keywords:',
                       help='List of IDs /keywords to be searched in NCBI databases')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        searchIds = []
        inputIds = self.getInputIds()

        for key in inputIds:
            searchIds.append(self._insertFunctionStep('searchStep', key, inputIds[key], prerequisites=[]))

        self._insertFunctionStep('createOutputStep', prerequisites=searchIds)

    def searchStep(self, key, maxEntries):
        dbName = self.getEnumText('dbType').lower()

        if self.searchMode.get() == IDS:
            ncbiIDs = [key]
        else:
            with Entrez.esearch(db=dbName, term=key, retmax=maxEntries, retmode='json') as handle:
                jDic = json.loads(handle.read())
                ncbiIDs = jDic['esearchresult']['idlist']

        if self.dbType.get() != 2:
            self.fetchSequences(ncbiIDs, dbName)
        else:
            self.fetchCompounds(ncbiIDs)

    def createOutputStep(self):
        if self.dbType.get() != 2:
            outputSet = SetOfSequences.create(self._getPath())
            outDir = self._getPath('sequences')
            for outFile in os.listdir(outDir):
                newSeq, outFile = Sequence(), os.path.abspath(os.path.join(outDir, outFile))
                newSeq.importFromFile(outFile, isAmino=self.dbType.get() == 0)
                outputSet.append(newSeq)
            self._defineOutputs(outputSequences=outputSet)

        else:
            outputSet = SetOfSmallMolecules(filename=self._getPath('outputSmallMolecules.sqlite'))
            outDir = self._getPath('compounds')
            for outFile in os.listdir(outDir):
                outMol = SmallMolecule(smallMolFilename=os.path.join(outDir, outFile))
                outMol.setMolName(os.path.splitext(os.path.basename(outFile))[0])
                outputSet.append(outMol)

            self._defineOutputs(outputSmallMolecules=outputSet)




    def fetchSequences(self, ncbiIDs, dbName):
        outDir = self._getPath('sequences')
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        for ncbiId in ncbiIDs:
            with Entrez.efetch(db=dbName, id=ncbiId, rettype='FASTA') as handle:
                with open(os.path.join(outDir, ncbiId+'.fa'), 'w') as f:
                    f.write(handle.read())

    def fetchCompounds(self, ncbiIDs):
        outDir = self._getPath('compounds')
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        for pID in ncbiIDs:
            outFile = os.path.abspath(os.path.join(outDir, '{}.sdf'.format(pID)))
            try:
                self.runJob('wget', '"{}" -O {}'.format(self.getPubChemURL(pID), outFile),
                            cwd=self._getTmpPath())
            except:
                try:
                    self.runJob('wget', '"{}" -O {}'.format(self.getPubChemURL(pID, dim=2), outFile),
                                cwd=self._getTmpPath())
                except:
                    print('Pubchem Compound with ID: {} could not be downloaded'.format(pID))


    def getPubChemURL(self, pID, dim=3):
        return "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{}/record/SDF/?record_type={}d&response_type=save&response_basename=Conformer{}D_CID_{}".\
            format(pID, dim, dim, pID)

    def getInputIds(self):
        ids = {}
        listIDs = self.listIDs.get()
        if self.searchMode.get() == IDS:
            for line in listIDs.split('\n'):
                if line.strip():
                    ids[json.loads(line.strip())['ID']] = 1

        elif self.searchMode.get() == KEYS:
            for line in listIDs.split('\n'):
                if line.strip():
                    jDic = json.loads(line.strip())
                    ids[jDic['ID']] = jDic['maxEntries']

        return ids

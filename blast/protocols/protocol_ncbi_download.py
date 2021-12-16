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

import os

from pwem.protocols import EMProtocol
from pwem.objects import Sequence, SetOfSequences
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pyworkflow.protocol.params import PointerParam, BooleanParam, StringParam, EnumParam, FileParam, STEPS_PARALLEL
from blast import Plugin
from pwchem.utils import guessIsAminoacids, parseFasta

class ProtChemNCBIDownload(EMProtocol):
    """Download the Fasta file(s) from NCBI databases"""
    _label = 'NCBI download'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('fromSet', BooleanParam, default=False,
                      label='Input from a set of databaseIds: ')
        form.addParam('inputIDSet', FileParam,
                       label='List of NCBI Ids:', condition='fromSet',
                       help="List of NCBI IDs for the query")
        form.addParam('inputID', StringParam,
                      label='NCBI Id:', allowsNull=False, condition='not fromSet',
                      help="NCBI ID for the query")

        form.addParam('dbType', EnumParam, default=0,
                      choices=['Protein', 'Nucleotide', 'SmallMolecule'], display=EnumParam.DISPLAY_HLIST,
                      label='Type of database to query on: ')
        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        searchIds = []
        if self.fromSet:
            inputIds = self.getInputIds()
            for ncbiID in inputIds:
                if self.dbType.get() != 2:
                    searchIds.append(self._insertFunctionStep('searchSequenceStep', ncbiID, prerequisites=[]))
                else:
                    searchIds.append(self._insertFunctionStep('searchCompoundStep', ncbiID, prerequisites=[]))
        else:
            if self.dbType.get() != 2:
                searchIds.append(self._insertFunctionStep('searchSequenceStep', self.inputID.get(), prerequisites=[]))
            else:
                searchIds.append(self._insertFunctionStep('searchCompoundStep', self.inputID.get(), prerequisites=[]))

        self._insertFunctionStep('createOutputStep', prerequisites=searchIds)

    def searchSequenceStep(self, ncbiID):
        dbName = self.getEnumText('dbType').lower()
        outDir = self._getPath('sequences')
        if not os.path.exists(outDir):
            os.mkdir(outDir)

        outFasta = os.path.abspath(('{}.fasta'.format(os.path.join(outDir, ncbiID))))
        fullCMD = '{} -db {} -id {} -format fasta > {}'.\
          format(Plugin.getEDirectProgram('efetch'), dbName, ncbiID, outFasta)

        Plugin.runEDirect(self, fullCMD, cwd=self._getPath())
        self.checkNotEmpty(outFasta, dbName, ncbiID)

    def searchCompoundStep(self, ncbiID):
        dbName = 'pccompound'
        outDir = self._getPath('compounds')
        if not os.path.exists(outDir):
            os.mkdir(outDir)

        outTxt = os.path.abspath(self._getTmpPath(ncbiID)) + '.txt'
        fullCMD = '{} -db {} -query {} | {} > {}'. \
            format(Plugin.getEDirectProgram('esearch'), dbName, ncbiID, Plugin.getEDirectProgram('efetch'), outTxt)

        Plugin.runEDirect(self, fullCMD, cwd=self._getPath())
        self.fetchFromPubChem(outTxt, ncbiID)

    def createOutputStep(self):
        if self.dbType.get() != 2:
            if self.fromSet:
                outputSet = SetOfSequences.create(self._getPath())
                outDir = self._getPath('sequences')
                fastaDic = {}
                for outFile in os.listdir(outDir):
                    fastaDic.update(parseFasta(os.path.join(outDir, outFile)))

                for seqId in fastaDic:
                    isAmino = guessIsAminoacids(fastaDic[seqId]['sequence'])
                    newSeq = Sequence(name=fastaDic[seqId]['name'], sequence=fastaDic[seqId]['sequence'],
                                      id=seqId, isAminoacids=isAmino)
                    outputSet.append(newSeq)
                self._defineOutputs(outputSequences=outputSet)
            else:
                outDir = self._getPath('sequences')
                for outFile in os.listdir(outDir):
                    #There should be just one
                    fastaDic = parseFasta(os.path.join(outDir, outFile), combined=False)
                    for seqId in fastaDic:
                        isAmino = guessIsAminoacids(fastaDic[seqId]['sequence'])
                        newSeq = Sequence(name=fastaDic[seqId]['name'], sequence=fastaDic[seqId]['sequence'],
                                          id=seqId, isAminoacids=isAmino)
                        self._defineOutputs(outputSequence=newSeq)

        else:
            outputSet = SetOfSmallMolecules(filename=self._getPath('outputSmallMolecules.sqlite'))
            outDir = self._getPath('compounds')
            for outFile in os.listdir(outDir):
                outMol = SmallMolecule(smallMolFilename=os.path.join(outDir, outFile))
                outputSet.append(outMol)

            self._defineOutputs(outputSmallMolecules=outputSet)

    def _validate(self):
        errors = []
        return errors

    def getDBName(self, dbText):
        return dbText.split('(')[-1].split(')')[0]

    def checkNotEmpty(self, file, dbName, ncbiID):
        delete = False
        with open(file) as f:
            if not f.read():
                delete = True
                print('\nNo sequence was found in database {} with query {}\n'.format(dbName, ncbiID))

        if delete:
            os.remove(file)

    def fetchFromPubChem(self, ncbiTxt, ncbiID):
        if os.path.getsize(ncbiTxt) > 0:
            with open(ncbiTxt) as fIn:
                for line in fIn:
                    if line.startswith('CID:'):
                        pID = line.strip().split()[1]

            outFile = os.path.abspath(self._getPath('compounds/{}.sdf'.format(pID)))
            try:
                self.runJob('wget', '"{}" -O {}'.format(self.getPubChemURL(pID), outFile),
                            cwd=self._getTmpPath())
            except:
                try:
                    self.runJob('wget', '"{}" -O {}'.format(self.getPubChemURL(pID, dim=2), outFile),
                                cwd=self._getTmpPath())
                except:
                    print('Pubchem Compound with ID: {} could not be downloaded'.format(pID))
        else:
            print('Your compound with ID: {} could not be found'.format(ncbiID))


    def getPubChemURL(self, pID, dim=3):
        return "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{}/record/SDF/?record_type={}d&response_type=save&response_basename=Conformer{}D_CID_{}".\
            format(pID, dim, dim, pID)

    def getInputIds(self):
        ids = []
        with open(self.inputIDSet.get()) as fIn:
            for line in fIn:
                ids.append(line.strip())
        return ids

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

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import *
from pyworkflow import BETA
from blast import Plugin
from pwem.objects import Sequence, SetOfSequences

from ..constants import *
from pwchem.utils import getSequenceFastaName

PROTEIN, NUCLEOTIDE = 0, 1

class ProtChemBLAST(EMProtocol):
    """Perform a BLAST search"""
    _label = 'BLAST search'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSequence', PointerParam, pointerClass='Sequence',
                      label='Input Sequence: ', allowsNull=False,
                      help="Sequence to be used as query")

        group.addParam('seqType', EnumParam, default=0,
                      choices=['Protein', 'Nucleotide'], display=EnumParam.DISPLAY_HLIST,
                      label='Type of sequence: ')

        group = form.addGroup('Database')
        group.addParam('localSearch', BooleanParam, default=False,
                      label='Local search: ')
        group.addParam('dbProtein', EnumParam, default=0,
                      choices=dbProtChoices, condition='not localSearch and seqType=={}'.format(PROTEIN),
                      label='Protein database to query on: ',
                      help='Nucleotide database to search on')
        group.addParam('dbNucleotide', EnumParam, default=0,
                      choices=dbNucChoices, condition='not localSearch and seqType=={}'.format(NUCLEOTIDE),
                      label='Nucleotide database to query on: ',
                      help='Nucleotide database to search on')

        group.addParam('dbName', EnumParam,
                       choices=Plugin.getLocalDatabases(),
                       label='Local database name: ', condition='localSearch',
                       help='Choose a database from those downloaded in {}'.format(Plugin.getDatabasesDir()))
        group.addParam('updateDB', BooleanParam, default=False, expertLevel=LEVEL_ADVANCED,
                       label='Update database: ', condition='localSearch',
                       help='In the case of being an NCBI database, update it before using it')

        group = form.addGroup('Program')
        group.addParam('blastProtein', EnumParam, default=0,
                      choices=self.getBLASTChoices(), condition='seqType=={}'.format(PROTEIN),
                      label='Protein BLAST type: ',
                      help='BLAST type to execute ')
        group.addParam('blastNucleotide', EnumParam, default=0,
                      choices=self.getBLASTChoices(protein=False), condition='seqType=={}'.format(NUCLEOTIDE),
                      label='Nucleotide BLAST type: ',
                      help='BLAST type to execute ')

        group.addParam('blastProteinProgram', EnumParam, default=0,
                      choices=self.getBLASTProgramChoices(),
                      condition='seqType=={} and blastProtein==0'.format(PROTEIN),
                      label='Protein BLAST program: ',
                      help='Protein BLAST program to execute:\n{}'.format(blastpProgramsHelp))
        group.addParam('blastNucleotideProgram', EnumParam, default=0,
                      choices=self.getBLASTProgramChoices(protein=False),
                      condition='seqType=={} and blastNucleotide==0'.format(NUCLEOTIDE),
                      label='Nucleotide BLAST program: ',
                      help='Nucleotide BLAST program to execute:\n{}'.format(blastnProgramsHelp))

        group.addParam('maxEntries', IntParam, default=20,
                       label='Maximum number of entries to keep: ',
                       help='Maximum number of entries to keep: -max_target_seqs')


        form.addSection(label='Parameters')
        form. addParam('labelAdvice', LabelParam,
                       label='BLAST default parameter will be used for each case if left empty.\n'
                             'Click to show default parameters for the defined input configuration',
                       help='For more information: '
                            'https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a')
        group = form.addGroup('General parameters')
        group.addParam('evalue', StringParam, default='0.05',
                      label='EValue for keep hits: ',
                      help='Expectation value (E) threshold for saving hits.\nIf empty, default will be used')
        group.addParam('word_size', StringParam, default='28',
                       label='Word size: ',
                       help='Word size for wordfinder algorithm (length of best perfect match).\n'
                            'If empty, default will be used')

        group = form.addGroup('Scoring parameters')
        #Matrix only not show for blastn
        group.addParam('matrix', EnumParam, default=4,
                       condition='not (seqType=={} and blastNucleotide==0)'.format(NUCLEOTIDE),
                       label='Scoring matrix: ', choices=matrixChoices,
                       help='Assigns a score for aligning pairs of residues, and determines overall alignment score')
        group.addParam('reward', StringParam, default='1',
                       label='Element match reward: ',
                       condition='seqType=={} and blastNucleotide==0'.format(NUCLEOTIDE),
                       help='Reward for a nucleotide / aminoacid match.\nIf empty, default will be used')
        group.addParam('penalty', StringParam, default='-2',
                       label='Element mismatch penalty: ',
                       condition='seqType=={} and blastNucleotide==0'.format(NUCLEOTIDE),
                       help='Penalty for a nucleotide / aminoacid mismatch.\nIf empty, default will be used')

        group.addParam('gapopen', StringParam, default='11',
                       label='Gap open cost: ',
                       condition='not (seqType=={} and blastNucleotide==2)'.format(PROTEIN),
                       help='Cost to open a gap.\nIf empty, default will be used')
        group.addParam('gapextend', StringParam, default='1',
                       label='Gap extend cost: ',
                       condition='not (seqType=={} and blastNucleotide==2)'.format(PROTEIN),
                       help='Cost to extend a gap.\nIf empty, default will be used')


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        searchIds = []
        searchIds.append(self._insertFunctionStep('BLASTSearchStep'))
        self._insertFunctionStep('createOutputStep')

    def BLASTSearchStep(self):
        if not self.localSearch.get():
            if self.seqType.get() == PROTEIN:
                dbName = self.getDBName(self.getEnumText('dbProtein'))
            else:
                dbName = self.getDBName(self.getEnumText('dbNucleotide'))
        else:
            dbName = self.getEnumText('dbName')

        outDir = self._getPath('sequences')
        if not os.path.exists(outDir):
            os.mkdir(outDir)

        inSeq = self.inputSequence.get()
        inFasta = os.path.abspath(self._getExtraPath(getSequenceFastaName(inSeq)))
        inSeq.exportToFile(inFasta)

        outFile = os.path.abspath(self._getPath(getSequenceFastaName(inSeq) + '.txt'))

        args = '-query {} -db {} -max_target_seqs {} -out {} -outfmt 4'.\
          format(inFasta, dbName, self.maxEntries.get(), outFile)
        if not self.localSearch.get():
            args += ' -remote'
        elif self.localSearch.get() and self.updateDB.get():
            upArgs = ' --decompress {} -passive'.format(dbName)
            Plugin.updateDatabase(self, upArgs)
        args += self.parseParameters()

        subprogram = self.getSelectedBLASTProgram()
        if subprogram in ['blastx', 'tblastn', 'tblastx']:
            program = subprogram
        elif subprogram in ['psi-blast', 'delta-blast']:
            program = subprogram.replace('-', '')
        elif subprogram in ['blastp', 'blastp-fast']:
            program = 'blastp'
            args += ' -task {}'.format(subprogram)
        elif subprogram in ['blastn', 'megablast', 'dc-megablast']:
            program = 'blastn'
            args += ' -task {}'.format(subprogram)

        Plugin.runBLAST(self, program, args, cwd=Plugin.getDatabasesDir())

    def createOutputStep(self):
        seqDic = self.parseBLASTOutput()
        outSeqs = SetOfSequences.create(self._getPath())

        #Adding input sequence
        inSeq = self.inputSequence.get()
        inSeq.evalue = Float(0.0)
        inSeq.score = Float(0.0)
        inSeq.firstPosition = Integer(1)
        outSeqs.append(inSeq)
        #Adding target sequences
        for seqId in seqDic:
            if seqId != getSequenceFastaName(inSeq):
                newSequence = seqDic[seqId]['sequence'].replace('-', '')
                isAmino = self.seqType.get() == 0
                newSeq = Sequence(name=seqId, sequence=newSequence, id=seqId, isAminoacids=isAmino,
                                  description=seqDic[seqId]['description'])
                newSeq.evalue = Float(seqDic[seqId]['evalue'])
                newSeq.score = Float(seqDic[seqId]['score'])
                newSeq.firstPosition = Integer((seqDic[seqId]['firstPosition']))
                outSeqs.append(newSeq)

        self._defineOutputs(outputSequences=outSeqs)


    def _validate(self):
        errors = []
        #Check numerical parameters
        for attr in self.getConditionalParameters():
            try:
                if getattr(self, attr).get() != '':
                    float(getattr(self, attr).get())
            except:
                errors.append('{} should be a number or empty'.format(attr))

        return errors

    def _warnings(self):
        warns = []
        # Check BLAST accepted parameter values
        if self.checkMatchMismatchType() == MATCH:
            # if blastn
            curMM = '{}/{}'.format(self.reward.get(), self.penalty.get())
            if curMM in blastnMM:
                gapPair = '{}/{}'.format(self.gapopen.get(), self.gapextend.get())
                if not gapPair in ALLOWED_BLAST_GAPS[curMM] and gapPair != '/':
                    warns.append('BLAST+ is somehow picky with the gap penalty values.\n'
                                  'Gap penalties {}/{} for match/mismatch values {} might yield an error.\n'
                                  'Documented options of gap penalties for these match/mismatch values are {}:'.
                                  format(self.gapopen.get(), self.gapextend.get(), curMM, ALLOWED_BLAST_GAPS[curMM]))
            else:
                warns.append('BLAST+ is somehow picky with the match / mismatch values.\n'
                            'Documented options of match/mismatch for blastn are: {}'.format(blastnMM))

        elif self.checkMatchMismatchType() == MATRIX:
            # if blastp, tblastn, blastx (not tblastx)
            curMatrix = self.getEnumText('matrix')
            gapPair = '{}/{}'.format(self.gapopen.get(), self.gapextend.get())
            if not gapPair in ALLOWED_BLAST_GAPS[curMatrix] and gapPair != '/':
                warns.append('BLAST+ is somehow picky with the gap penalty values.\n'
                             'Gap penalties {}/{} for {} matrix might yield an error.\n'
                             'Documented options of gap penalties for this matrix are {}:'.
                             format(self.gapopen.get(), self.gapextend.get(), curMatrix, ALLOWED_BLAST_GAPS[curMatrix]))

        if warns != []:
            warns.append('\n\nDo you want to keep trying to BLAST with the specified parameters?')
        return warns

#################### UTILS ##################

    def getDBChoices(self):
        if self.seqType.get() == PROTEIN:
            return dbProtChoices
        elif self.seqType.get() == NUCLEOTIDE:
            return dbNucChoices

    def getBLASTChoices(self, protein=True):
        if protein:
            return ['blastp', 'tblastn']
        else:
            return ['blastn', 'blastx', 'tblastx']

    def getBLASTProgramChoices(self, protein=True):
        if protein:
            return ['blastp', 'blastp-fast', 'psi-blast', 'delta-blast']
        else:
            return ['blastn', 'megablast', 'dc-megablast']

    def getDBName(self, dbText):
        return dbText.split('(')[-1].split(')')[0]

    def checkMatchMismatchType(self):
        if self.seqType.get() == NUCLEOTIDE and self.blastNucleotide.get() == 0:
            #Blastn: match/mismatch penalties and gap penalties
            return MATCH
        elif self.seqType.get() == PROTEIN or (self.seqType.get() == NUCLEOTIDE and self.blastNucleotide.get() != 2):
            #blastp, tblastn, blastx: penalty matrix and gap penalties
            return MATRIX
        else:
            #tblastx: penalty matrix, no gap penalty
            return NOGAP
        
    def getSelectedBLASTProgram(self):
        if self.seqType.get() == PROTEIN:
            if self.blastProtein.get() == 0:
                #BLASTP subprogram
                return self.getEnumText('blastProteinProgram')
            else:
                #tblastn
                return self.getEnumText('blastProtein')
        else:
            if self.blastNucleotide.get() == 0:
                #BLASTN subprogram
                return self.getEnumText('blastNucleotideProgram')
            else:
                #'blastx', 'tblastx'
                return self.getEnumText('blastNucleotide')

    #PARAMETERS PARSING
    def parseParameters(self):
        parArgs = ''
        for parName in self.getConditionalParameters():
            if getattr(self, parName).get() != '':
                parArgs += ' -{} {}'.format(parName, getattr(self, parName).get())

        if self.checkMatchMismatchType() != MATCH:
            parArgs += ' -matrix {}'.format(self.getEnumText('matrix'))

        return parArgs

    def getConditionalParameters(self):
        if self.checkMatchMismatchType() == MATCH:
            return ['evalue', 'word_size', 'reward', 'penalty', 'gapopen', 'gapextend']
        elif self.checkMatchMismatchType() == MATRIX:
            return ['evalue', 'word_size', 'gapopen', 'gapextend']
        else:
            return ['evalue', 'word_size']


    def parseBLASTOutput(self):
        def goToNextLine(fIn, read=None):
            line = ''
            while line.strip() == '':
                line = fIn.readline()
            if read is not None:
                read += 1
            return line, read

        inSeq = self.inputSequence.get()
        outFile = os.path.abspath(self._getPath(getSequenceFastaName(inSeq) + '.txt'))
        seqDic, read = {}, 0
        with open(outFile) as fIn:
            for line in fIn:
                if read == 0 and line.startswith('Sequences producing significant alignments'):

                    read, line = 1, fIn.readline()
                    while line.strip() == '':
                        line = fIn.readline()

                if read == 1:
                    if line.strip() != '':
                        line = line.strip().split()
                        seqDic[line[0]] = {'description': ' '.join(line[1:-2]), 'score': line[-2], 'evalue': line[-1],
                                           'sequence': ''}
                    else:
                        line, read = goToNextLine(fIn, read)

                if read == 2:
                    if line.strip() != '':
                        line = line.strip().split()
                        if line[0] in seqDic:
                            seqDic[line[0]]['sequence'] += line[2]
                            seqDic[line[0]]['firstPosition'] = line[1]
                    else:
                        line, read = goToNextLine(fIn, read)

                if read == 3:
                    if line.strip() != '':
                        line = line.strip().split()
                        if line[0] in seqDic:
                            seqDic[line[0]]['sequence'] += line[2]
        return seqDic













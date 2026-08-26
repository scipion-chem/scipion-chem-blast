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

from pyworkflow.protocol.params import *
from pyworkflow import BETA
from pwem.objects import Sequence, SetOfSequences
from pwem.protocols import EMProtocol

from pwchem.utils import getSequenceFastaName
from pwchem.utils.utilsFasta import fastFastaExport

from blast import Plugin
from blast.constants import *


PROTEIN, NUCLEOTIDE = 0, 1
TIERED_EVALUE_CONDITION = 'multipleQueries and autoTieredEvalue'

class ProtChemBLAST(EMProtocol):
    """Perform a BLAST search.

    AI Generated: batch mode support (multipleQueries) was added to accept a whole
    SetOfSequences/SetOfSequenceROIs as query instead of a single Sequence, running one BLAST
    call per input set (or per length tier if "Auto length-tiered E-value" is enabled), useful
    for any workflow needing to BLAST many short queries against the same database without one
    protocol run per query.
    """
    _label = 'BLAST search'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('multipleQueries', BooleanParam, default=False,
                      label='Batch multiple queries: ',
                      help='Instead of a single Sequence, BLAST a whole SetOfSequences or '
                           'SetOfSequenceROIs in one batch. All queries are grouped into one BLAST '
                           'call, or several length-based tiers if "Auto length-tiered E-value" below '
                           'is enabled.')
        group.addParam('inputSequence', PointerParam, pointerClass='Sequence',
                      label='Input Sequence: ', allowsNull=False, condition='not multipleQueries',
                      help="Sequence to be used as query")
        group.addParam('inputSequences', PointerParam, pointerClass='SetOfSequences,SetOfSequenceROIs',
                      label='Input Sequences: ', allowsNull=True, condition='multipleQueries',
                      help="Set of sequences/ROIs to be used as query")

        group.addParam('seqType', EnumParam, default=1,
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
                       help='Maximum number of entries to keep: -max_target_seqs. '
                            'Undefined with <0')


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

        tGroup = form.addGroup('Batch length tiers', condition='multipleQueries')
        tGroup.addParam('autoTieredEvalue', BooleanParam, default=False, condition='multipleQueries',
                       label='Auto length-tiered E-value: ',
                       help='Group queries into short/medium/long length tiers and BLAST each tier with '
                            'its own E-value instead of the single one above. Standard BLAST statistics '
                            'penalize short queries: a real short hit against a large database can be '
                            'discarded as "not significant" under a strict default E-value, so a laxer '
                            'threshold for short queries avoids losing real short matches.')
        tGroup.addParam('shortMaxLen', IntParam, default=30, condition=TIERED_EVALUE_CONDITION,
                       expertLevel=LEVEL_ADVANCED, label='Short query max length: ',
                       help='Query length (aa/nt) up to which the "short" E-value tier is used.')
        tGroup.addParam('mediumMaxLen', IntParam, default=100, condition=TIERED_EVALUE_CONDITION,
                       expertLevel=LEVEL_ADVANCED, label='Medium query max length: ',
                       help='Query length (aa/nt) up to which the "medium" E-value tier is used (above it, '
                            '"long" is used).')
        tGroup.addParam('evalueShort', StringParam, default='50', condition=TIERED_EVALUE_CONDITION,
                       expertLevel=LEVEL_ADVANCED, label='E-value (short): ')
        tGroup.addParam('evalueMedium', StringParam, default='0.1', condition=TIERED_EVALUE_CONDITION,
                       expertLevel=LEVEL_ADVANCED, label='E-value (medium): ')
        tGroup.addParam('evalueLong', StringParam, default='0.05', condition=TIERED_EVALUE_CONDITION,
                       expertLevel=LEVEL_ADVANCED, label='E-value (long): ')

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
        dbName = self.resolveDBName()

        outDir = self._getPath('sequences')
        if not os.path.exists(outDir):
            os.mkdir(outDir)

        if self.localSearch.get() and self.updateDB.get():
            upArgs = ' --decompress {} -passive'.format(dbName)
            Plugin.updateDatabase(self, upArgs)

        if not self.multipleQueries.get():
            inSeq = self.inputSequence.get()
            inFasta = os.path.abspath(self._getExtraPath(getSequenceFastaName(inSeq) + '.fasta'))
            inSeq.exportToFile(inFasta)
            outFile = os.path.abspath(self._getPath(getSequenceFastaName(inSeq) + '.txt'))
            self.runBlastCall(inFasta, dbName, outFile)
        else:
            for tierIdx, (evalueOverride, tierQueries) in enumerate(self.getBatchTiers()):
                if not tierQueries:
                    continue
                inFasta = self.getBatchTierFastaFile(tierIdx)
                with open(inFasta, 'w') as fh:
                    for qId, qSeq in tierQueries:
                        fh.write('>{}\n{}\n'.format(qId, qSeq))
                outFile = self.getBatchTierOutFile(tierIdx)
                self.runBlastCall(inFasta, dbName, outFile, evalueOverride=evalueOverride)

    def runBlastCall(self, inFasta, dbName, outFile, evalueOverride=None):
        args = '-query {} -db {} -out {} -outfmt 15'.format(inFasta, dbName, outFile)
        if self.maxEntries.get() > 0:
            args += ' -max_target_seqs {}'.format(self.maxEntries.get())
        if not self.localSearch.get():
            args += ' -remote'
        args += self.parseParameters(evalueOverride=evalueOverride)

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
        outSeqs = SetOfSequences.create(self._getPath())

        if not self.multipleQueries.get():
            self._appendSingleQueryOutput(outSeqs)
        else:
            self._appendBatchOutput(outSeqs)

        outPath = self._getExtraPath('viewSequences.fasta')
        fastFastaExport(outSeqs, outPath)
        # outSeqs.exportToFile(outPath)

        self._defineOutputs(outputSequences=outSeqs)

    def _appendSingleQueryOutput(self, outSeqs):
        inSeq = self.inputSequence.get()
        outFile = os.path.abspath(self._getPath(getSequenceFastaName(inSeq) + '.txt'))
        seqDic = self.parseBLASTOutput(outFile)
        isAmino = self.seqType.get() == 0

        #Adding target sequences
        for seqId in seqDic:
            if seqId != 'Query_1':
                newSequence = seqDic[seqId]['sequence']
                newSeq = Sequence(name=seqId, sequence=newSequence, id=seqId, isAminoacids=isAmino,
                                  description=seqDic[seqId]['description'])
                newSeq.evalue = String(str(seqDic[seqId]['evalue']))
                newSeq.score = Float(seqDic[seqId]['score'])
                outSeqs.append(newSeq)
            else:
                inSeq.setSequence(seqDic[seqId]['sequence'])
                inSeq.evalue = Float(0.0)
                inSeq.score = Float(0.0)
                outSeqs.append(inSeq)

    def _appendBatchOutput(self, outSeqs):
        isAmino = self.seqType.get() == 0
        for tierIdx, (_, tierQueries) in enumerate(self.getBatchTiers()):
            if not tierQueries:
                continue
            outFile = self.getBatchTierOutFile(tierIdx)
            seqDic = self.parseBLASTOutput(outFile)
            for seqId, info in seqDic.items():
                newSeq = Sequence(name=seqId, sequence=info['sequence'], id=seqId, isAminoacids=isAmino,
                                  description=info['description'])
                newSeq.evalue = String(str(info['evalue']))
                newSeq.score = Float(info['score'])
                newSeq.queryId = String(info['query_title'])
                outSeqs.append(newSeq)


    def _validate(self):
        errors = self._validateConditionalParameters()

        if self.seqType.get() == 0 and int(self.word_size.get()) >= 8:
            errors.append('Word size must be < 8 when using blastp. Check the specified parameters.')

        if not self.multipleQueries.get() and self.inputSequence.get() is None:
            errors.append('You must provide an Input Sequence.')
        if self.multipleQueries.get() and self.inputSequences.get() is None:
            errors.append('You must provide an Input Sequences set.')

        if self.multipleQueries.get() and self.autoTieredEvalue.get():
            errors += self._validateTieredEvalue()

        return errors

    def _validateConditionalParameters(self):
        errors = []
        #Check numerical parameters
        for attr in self.getConditionalParameters():
            try:
                if getattr(self, attr).get() != '':
                    float(getattr(self, attr).get())
            except:
                errors.append('{} should be a number or empty'.format(attr))
        return errors

    def _validateTieredEvalue(self):
        errors = []
        for attr in ('evalueShort', 'evalueMedium', 'evalueLong'):
            try:
                float(getattr(self, attr).get())
            except (TypeError, ValueError):
                errors.append('{} should be a number'.format(attr))
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

    def resolveDBName(self):
        if not self.localSearch.get():
            if self.seqType.get() == PROTEIN:
                return self.getDBName(self.getEnumText('dbProtein'))
            else:
                return self.getDBName(self.getEnumText('dbNucleotide'))
        else:
            return self.getEnumText('dbName')

    def getBatchQueries(self):
        '''Return a list of (queryId, sequence) tuples from inputSequences, which can be a
        SetOfSequences or a SetOfSequenceROIs.'''
        queries = []
        for item in self.inputSequences.get():
            item = item.clone()
            if hasattr(item, 'getROISequence'):
                queries.append((item.getROIId(), item.getROISequence()))
            else:
                queries.append((item.getSeqName() or item.getId(), item.getSequence()))
        return queries

    def getBatchTiers(self):
        '''Group batch queries into (evalueOverride, queryList) tuples. With auto length-tiered
        E-value disabled, returns a single tier using the general "evalue" parameter (evalueOverride
        left as None). Otherwise splits queries into short/medium/long length tiers, each with its
        own E-value.'''
        queries = self.getBatchQueries()
        if not self.autoTieredEvalue.get():
            return [(None, queries)]

        shortQ, mediumQ, longQ = [], [], []
        shortMax, mediumMax = self.shortMaxLen.get(), self.mediumMaxLen.get()
        for qId, qSeq in queries:
            n = len(qSeq)
            if n <= shortMax:
                shortQ.append((qId, qSeq))
            elif n <= mediumMax:
                mediumQ.append((qId, qSeq))
            else:
                longQ.append((qId, qSeq))
        return [
            (self.evalueShort.get(), shortQ),
            (self.evalueMedium.get(), mediumQ),
            (self.evalueLong.get(), longQ),
        ]

    def getBatchTierFastaFile(self, tierIdx):
        return os.path.abspath(self._getExtraPath('batch_tier{}.fasta'.format(tierIdx)))

    def getBatchTierOutFile(self, tierIdx):
        return os.path.abspath(self._getPath('batch_tier{}.txt'.format(tierIdx)))

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
    def parseParameters(self, evalueOverride=None):
        parArgs = ''
        for parName in self.getConditionalParameters():
            val = getattr(self, parName).get()
            if parName == 'evalue' and evalueOverride is not None:
                val = evalueOverride
            if val != '':
                parArgs += ' -{} {}'.format(parName, val)

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


    def parseBLASTOutput(self, outFile):
      '''Parse a BLAST JSON (-outfmt 15) output file. Iterates every query block in
      BlastOutput2 (a multi-record query FASTA produces one block per query), not just the
      first one, so this also supports batch/multi-query files.'''
      with open(outFile) as f:
        js = json.loads(f.read())

      blocks = js['BlastOutput2']
      multiQuery = len(blocks) > 1
      seqDic = {}
      for block in blocks:
        results = block['report']['results']['search']
        queryTitle = results.get('query_title') or results.get('query_id', '')
        hits = results['hits']
        for h in hits:
          accession, hId = h['description'][0]['accession'], h['description'][0]['id']
          for hs in h['hsps']:
            hf, qf, ql, ev, sc = hs['hit_from'], hs['query_from'], hs['query_to'], hs['evalue'], hs['score']
            hseq = hs['hseq']
            acc = f'{accession}_{hf}'
            seqId = f'{queryTitle}::{acc}' if multiQuery else acc

            seqDic[seqId] = {'score': sc, 'evalue': ev, 'sequence': hseq.replace('-', ''),
                              'description': hId, 'query_title': queryTitle}

      return seqDic













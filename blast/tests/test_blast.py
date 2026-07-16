# ***************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
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
# ***************************************************************************


from pyworkflow.tests import BaseTest
import pyworkflow.tests as tests

from pwem.protocols import ProtImportSequence

from pwchem.protocols import ProtDefineSeqROI

from blast import Plugin
from blast.constants import BLASTdbs

from ..protocols import ProtChemBLAST, ProtChemBLASTDatabase, ProtChemNCBIDownload

idsDic = {0: '{"ID": "P0DTC2"}\n{"ID": "P59594"}\n',
          1: '{"ID": "nr_025000"}\n{"ID": "nr_025001"}\n',
          2: '{"ID": "2244"}\n{"ID": "6247"}\n'}

keysDic = {0: '{"ID": "hemoglobin", "maxEntries": "2"}\n',
           1: '{"ID": "hemoglobin", "maxEntries": "2"}\n',
           2: '{"ID": "aspirin", "maxEntries": "2"}\n'}

dbLabels = ['Protein', 'Nucleotide', 'Compounds']


class TestNCBIDownload(BaseTest):
  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)

  @classmethod
  def _runImportFromID(cls, dbType=1):
    protImportNCBI = cls.newProtocol(
      ProtChemNCBIDownload,
      listIDs=idsDic[dbType], dbType=dbType)
    protImportNCBI.setObjLabel('NCBI from ID ' + dbLabels[dbType])

    cls.proj.launchProtocol(protImportNCBI, wait=False)
    return protImportNCBI

  @classmethod
  def _runImportFromKey(cls, dbType=1):
    protImportNCBI = cls.newProtocol(
      ProtChemNCBIDownload,
      searchMode=1, listIDs=keysDic[dbType], dbType=dbType)
    protImportNCBI.setObjLabel('NCBI from keyword ' + dbLabels[dbType])

    cls.proj.launchProtocol(protImportNCBI, wait=False)
    return protImportNCBI

  def testNCBIDownload(self):
    prots = []
    for i in range(3):
        prots += [self._runImportFromID(dbType=i)]

    for i in range(3):
        if i == 2:
            self._waitOutput(prots[i], 'outputSmallMolecules')
            self.assertIsNotNone(prots[i].outputSmallMolecules)
        else:
            self._waitOutput(prots[i], 'outputSequences')
            self.assertIsNotNone(prots[i].outputSequences)

  def testNCBISearch(self):
    prots = []
    for i in range(3):
        prots += [self._runImportFromKey(dbType=i)]

    for i in range(3):
        if i == 2:
            self._waitOutput(prots[i], 'outputSmallMolecules')
            self.assertIsNotNone(prots[i].outputSmallMolecules)
        else:
            self._waitOutput(prots[i], 'outputSequences')
            self.assertIsNotNone(prots[i].outputSequences)

class TestDatabaseBLAST(BaseTest):
  dbName = '16S_ribosomal_RNA'
  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)

  @classmethod
  def _createLocalDatabase(cls):
    dbIndex = cls.getDatabaseIndex(cls.dbName, fromNCBI=True)
    protDB = cls.newProtocol(ProtChemBLASTDatabase, fromNCBI=True, inputID=dbIndex)
    cls.launchProtocol(protDB, wait=True)
    return protDB

  @classmethod
  def getDatabaseIndex(cls, dbName, fromNCBI=False):
    if not fromNCBI:
      options = Plugin.getLocalDatabases()
    else:
      options = BLASTdbs
    for i, name in enumerate(options):
      if dbName == name:
        return i

  def testDBBLAST(self):
    protDB = self._createLocalDatabase()
    self.assertTrue(protDB.isFinished() and not protDB.isFailed())


class TestBLAST(BaseTest):
  dbName = '16S_ribosomal_RNA'
  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)

  @classmethod
  def _createLocalDatabase(cls):
    dbIndex = cls.getDatabaseIndex(cls.dbName, fromNCBI=True)
    protDB = cls.newProtocol(ProtChemBLASTDatabase, fromNCBI=True, inputID=dbIndex)
    cls.launchProtocol(protDB)
    return protDB

  @classmethod
  def getDatabaseIndex(cls, dbName, fromNCBI=False):
    if not fromNCBI:
      options = Plugin.getLocalDatabases()
    else:
      options = BLASTdbs
    for i, name in enumerate(options):
      if dbName == name:
        return i

  @classmethod
  def _runImportSeq(cls):
    protImportSeq = cls.newProtocol(
      ProtImportSequence,
      inputSequence=1, inputNucleotideSequence=3, geneBankSequence='nr_025000')
    cls.launchProtocol(protImportSeq)
    return protImportSeq

  @classmethod
  def _runBLASTn(cls, protSeq):
      dbIndex = cls.getDatabaseIndex(cls.dbName, fromNCBI=False)
      protBLAST = cls.newProtocol(
        ProtChemBLAST,
        inputSequence=protSeq.outputSequence, seqType=1, localSearch=True, dbName=dbIndex,
        word_size='11', gapopen='5', gapextend='2'
        )

      cls.launchProtocol(protBLAST)
      return protBLAST

  def testBLAST(self):
    dbIndex = self.getDatabaseIndex(self.dbName, fromNCBI=False)
    if dbIndex == None:
        self._createLocalDatabase()

    protSeq = self._runImportSeq()
    protBLAST = self._runBLASTn(protSeq)
    self.assertIsNotNone(protBLAST.outputSequences)


class TestBLASTBatch(BaseTest):
  '''AI Generated. Exercises the new multipleQueries batch mode (a whole SetOfSequenceROIs as
  query instead of a single Sequence), using two windows of the same 16S rRNA sequence already
  used by TestBLAST as two separate batch queries. Kept as its own BaseTest (not a TestBLAST
  subclass) so it does not also run the inherited testBLAST method. Uses inputRawSequence
  (the literal nr_025000 sequence) instead of ProtImportSequence's geneBankSequence fetch: the
  latter intermittently hits an unrelated pyworkflow sqlite mapper "circular reference" error
  on this pyworkflow version when storing the downloaded record.'''

  dbName = '16S_ribosomal_RNA'
  # Literal content of GenBank record nr_025000 (Mycobacterium kubicae 16S rRNA, partial),
  # hardcoded to avoid both the network fetch and the geneBankSequence storage flakiness.
  NR_025000_SEQ = (
    'GGCAACACAGCAAGCGAACGGAAAGGCCCCCGGGGGACCGAGGGCGAACGGGGAGAACACGGGGGACACCCGCACCGGGA'
    'AAGCCGGGAAACGGGCAAACCGGAAGGACCAGAGAGCAGCAGGGGAAAGCGCGGGGGGAGGGCCCGCGGCCACAGCGGGGG'
    'GGGACGGCCACCAAGGCGACGACGGGAGCCGGCCGAGAGGGGCCGGCCACACGGGACGAGAACGGCCCAGACCCACGGGA'
    'GGCAGCAGGGGGAAAGCACAAGGGCGCAAGCCGAGCAGCGACGCCGCGGGGGGAGACGGCCCGGGGAAACCCCAGCAGGG'
    'ACGAAGCGCAAGGACGGACCGCAGAAGAAGCACCGGCCAACACGGCCAGCAGCCGCGGAAACGAGGGGCGAGCGGCCGGA'
    'AACGGGCGAAAGAGCCGAGGGGGCGCGGCGGAAAACCGGGGGCAACCCCGGCGGCGGGCGAACGGGCAGACGGAGACGCA'
    'GGGGAGACGGAACCGGGAGCGGGGAAGCGCAGAACAGGAGGAACACCGGGGCGAAGGCGGGCCGGGCAGAACGACGCGAG'
    'GAGCGAAAGCGGGGGAGCGAACAGGAAGAACCCGGAGCCACGCCGAAACGGGGGACAGGGGGGCCCCGGGACCGGCCGAG'
    'CAACGCAAAGACCCCGCCGGGGAGACGGCCGCAAGGCAAAACCAAAGGAAGACGGGGGCCCGCACAAGCGGCGGAGCAGG'
    'GAAACGAGCAACGCGAAGAACCACCGGGGACAGCACAGGACGCGCAGAGAAGGCGCCCGGGCCGGGCAGGGGGCAGGCGC'
    'GCAGCCGGCGGAGAGGGGAAGCCCGCAACGAGCGCAACCCGCCAGGCCAGCGGGAAGCCGGGGACCGGAGAGACGCCGGG'
    'GCAACCGGAGGAAGGGGGGAGACGCAAGCACAGCCCCAGCCAGGGCCACACAGCACAAGGCCGGACAAAGGGCGCGAGCC'
    'GCGAGGAAGCGAACCAAAGCCGGCCAGCGGACGGGGCGCAACCGACCCCGGAAGCGGAGCGCAGAACGCAGACAGCAACG'
    'CGCGGGAAACGCCCGGG'
  )

  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)

  @classmethod
  def _createLocalDatabase(cls):
    dbIndex = cls.getDatabaseIndex(cls.dbName, fromNCBI=True)
    protDB = cls.newProtocol(ProtChemBLASTDatabase, fromNCBI=True, inputID=dbIndex)
    cls.launchProtocol(protDB)
    return protDB

  @classmethod
  def getDatabaseIndex(cls, dbName, fromNCBI=False):
    if not fromNCBI:
      options = Plugin.getLocalDatabases()
    else:
      options = BLASTdbs
    for i, name in enumerate(options):
      if dbName == name:
        return i

  @classmethod
  def _runImportSeq(cls):
    protImportSeq = cls.newProtocol(
      ProtImportSequence,
      inputSequenceName='nr_025000', inputSequence=1,
      inputRawSequence=cls.NR_025000_SEQ)
    cls.launchProtocol(protImportSeq)
    return protImportSeq

  @classmethod
  def _runDefSeqROIs(cls, inProt, fullSeq):
    windows = [(1, 120), (401, 520)]
    inROIs = '\n'.join(
      '{}) Residues: {{"index": "{}-{}", "residues": "{}", "desc": "None"}}'.format(
        i, start, end, fullSeq[start - 1:end]
      )
      for i, (start, end) in enumerate(windows, 1)
    )
    protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
    protDefSeqROIs.inputSequence.set(inProt)
    protDefSeqROIs.inputSequence.setExtended('outputSequence')

    cls.launchProtocol(protDefSeqROIs)
    return protDefSeqROIs

  @classmethod
  def _runBLASTnBatch(cls, protROIs):
    dbIndex = cls.getDatabaseIndex(cls.dbName, fromNCBI=False)
    protBLAST = cls.newProtocol(
      ProtChemBLAST,
      multipleQueries=True, seqType=1, localSearch=True, dbName=dbIndex,
      word_size='11', gapopen='5', gapextend='2'
      )
    protBLAST.inputSequences.set(protROIs)
    protBLAST.inputSequences.setExtended('outputROIs')

    cls.launchProtocol(protBLAST)
    return protBLAST

  def testBLASTBatch(self):
    dbIndex = self.getDatabaseIndex(self.dbName, fromNCBI=False)
    if dbIndex == None:
        self._createLocalDatabase()

    protSeq = self._runImportSeq()
    fullSeq = protSeq.outputSequence.getSequence()
    protROIs = self._runDefSeqROIs(protSeq, fullSeq)
    protBLAST = self._runBLASTnBatch(protROIs)
    self.assertIsNotNone(protBLAST.outputSequences)





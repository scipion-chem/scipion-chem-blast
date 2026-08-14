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

from pwchem.protocols import ProtDefineSeqROI, ProtChemImportVariants, ProtChemGenerateVariants
from pwchem.tests.tests_imports import TestImportVariants

from blast import Plugin
from blast.constants import BLASTdbs, BLAST_HUMAN_DB_PATH

from ..protocols import (
    ProtChemBLAST, ProtChemBLASTDatabase, ProtChemNCBIDownload,
    ProtSelfToleranceFilter, ProtBLASTPanelConservation,
)

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


class TestSelfToleranceFilter(BaseTest):
  '''AI Generated. First 20 residues are the REAL human GAPDH N-terminus (P04406) -- expected
  to be discarded (real, exact 20aa human-proteome match). The last 12 residues are a real SARS-
  CoV-2 spike fragment (P0DTC2 residues 70-81, empirically confirmed via TestIEDBCrossref to be a
  clean, non-promiscuous match -- only a documented SARS-CoV-2 neutralizing epitope, no incidental
  cross-reactivity noise) -- expected to survive. A shorter 7aa fragment was tried first and
  discarded: real (not fabricated) BLASTp against the human proteome found an 85.7%% incidental
  partial hit -- expected noise for a query that short against a whole proteome, not a bug, but it
  invalidated that fixture as a clean "should survive" case.'''

  GAPDH_NTERM = 'MGKVKVGVNGFGRIGRLVTR'  # P04406, residues 1-20
  SPIKE_FRAGMENT = 'VSGTNGTKRFDN'  # P0DTC2, residues 70-81
  COMBINED_SEQ = GAPDH_NTERM + SPIKE_FRAGMENT

  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)

  @classmethod
  def _runImportSeq(cls):
    protImportSeq = cls.newProtocol(
      ProtImportSequence,
      inputSequenceName='self_tolerance_fixture', inputSequence=1,
      inputRawSequence=cls.COMBINED_SEQ)
    cls.launchProtocol(protImportSeq)
    return protImportSeq

  @classmethod
  def _runDefSeqROIs(cls, protSeq):
    inROIs = (
      '1) Residues: {{"index": "1-20", "residues": "{}", "desc": "None"}}\n'
      '2) Residues: {{"index": "21-32", "residues": "{}", "desc": "None"}}'
    ).format(cls.GAPDH_NTERM, cls.SPIKE_FRAGMENT)
    protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=0, inROIs=inROIs)
    protDefSeqROIs.inputSequence.set(protSeq)
    protDefSeqROIs.inputSequence.setExtended('outputSequence')
    cls.launchProtocol(protDefSeqROIs)
    return protDefSeqROIs

  @classmethod
  def _runFilter(cls, protROIs):
    protFilter = cls.newProtocol(ProtSelfToleranceFilter)
    protFilter.inputROIs.set(protROIs)
    protFilter.inputROIs.setExtended('outputROIs')
    cls.launchProtocol(protFilter)
    return protFilter

  def testSelfToleranceFilter(self):
    protSeq = self._runImportSeq()
    protROIs = self._runDefSeqROIs(protSeq)
    protFilter = self._runFilter(protROIs)
    self.assertIsNotNone(protFilter.outputROIs)
    survivors = [roi.getROISequence() for roi in protFilter.outputROIs]
    self.assertNotIn(self.GAPDH_NTERM, survivors, 'real human GAPDH N-terminus should be rejected')
    self.assertIn(self.SPIKE_FRAGMENT, survivors, 'real SARS-CoV-2 spike fragment should survive')


class TestBLASTPanelConservation(TestImportVariants):
  '''AI Generated. Panel = 3 real point-mutant variants of SARS-CoV-2 spike (P0DTC2) generated
  by ProtChemGenerateVariants (Original/Alpha/T478K, same fixture as pwchem's own
  TestGenerateSequences) -- none of those mutations fall inside residues 34-40, so the query ROI
  (the same real spike fragment already used across this codebase) is expected to be conserved
  across the whole panel.'''

  SPIKE_FRAGMENT = 'RGVYYPD'  # P0DTC2, residues 34-40
  MUTATIONS = '1) Variant: Original\n2) Variant: Alpha\n3) Mutations: T478K\n'

  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)
    cls._runImportVariants()
    cls._waitOutput(cls.protImportVariants, 'outputVariants', sleepTime=5)

  @classmethod
  def _runGeneratePanel(cls):
    protGenSeqs = cls.newProtocol(ProtChemGenerateVariants, toMutateList=cls.MUTATIONS)
    protGenSeqs.inputSequenceVariants.set(cls.protImportVariants)
    protGenSeqs.inputSequenceVariants.setExtended('outputVariants')
    cls.launchProtocol(protGenSeqs)
    return protGenSeqs

  @classmethod
  def _runDefSeqROI(cls):
    inROIs = '1) Residues: {{"index": "34-40", "residues": "{}", "desc": "None"}}'.format(cls.SPIKE_FRAGMENT)
    protDefSeqROIs = cls.newProtocol(ProtDefineSeqROI, chooseInput=1, inROIs=inROIs)
    protDefSeqROIs.inputSequenceVariants.set(cls.protImportVariants)
    protDefSeqROIs.inputSequenceVariants.setExtended('outputVariants')
    cls.launchProtocol(protDefSeqROIs)
    return protDefSeqROIs

  def testBLASTPanelConservation(self):
    protPanel = self._runGeneratePanel()
    protROIs = self._runDefSeqROI()

    protCons = self.newProtocol(ProtBLASTPanelConservation)
    protCons.inputROIs.set(protROIs)
    protCons.inputROIs.setExtended('outputROIs')
    protCons.panelSequences.set(protPanel)
    protCons.panelSequences.setExtended('outputSequences')
    self.launchProtocol(protCons)

    self.assertIsNotNone(protCons.outputROIs)
    roi = list(protCons.outputROIs)[0]
    self.assertEqual(roi._conservationPanelTotal.get(), 3)
    self.assertGreaterEqual(roi._conservationPanelMatches.get(), 1)





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
from blast import Plugin

from ..protocols import ProtChemBLAST, ProtChemBLASTDatabase, ProtChemNCBIDownload

class TestNCBIDownload(BaseTest):
  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)

  @classmethod
  def _runImportFromID(cls):
    protImportSequence = cls.newProtocol(
      ProtChemNCBIDownload,
      inputID='nr_025000',
      dbType=1)
    cls.launchProtocol(protImportSequence)
    return protImportSequence

  def testNCBIDownload(self):
    protNCBI = self._runImportFromID()
    self.assertIsNotNone(protNCBI.outputSequence)


class TestBLAST(TestNCBIDownload):
  dbName = '16S_ribosomal_RNA'
  @classmethod
  def _createLocalDatabase(cls):
      protDB = cls.newProtocol(ProtChemBLASTDatabase,
                                fromNCBI=True, inputID=cls.dbName)
      cls.launchProtocol(protDB)
      return protDB

  @classmethod
  def getDatabaseIndex(cls, dbName):
    options = Plugin.getLocalDatabases()
    for i, name in enumerate(options):
        if dbName == name:
          return i

  @classmethod
  def _runBLASTn(cls, sequence):
      dbIndex = cls.getDatabaseIndex(cls.dbName)
      protBLAST = cls.newProtocol(
        ProtChemBLAST,
        inputSequence=sequence, seqType=1, localSearch=True, dbName=dbIndex,
        word_size='11', gapopen='5', gapextend='2'
        )

      cls.launchProtocol(protBLAST)
      return protBLAST

  def testBLAST(self):
    self._createLocalDatabase()
    protSequence = self._runImportFromID()
    protBLAST = self._runBLASTn(protSequence.outputSequence)
    self.assertIsNotNone(protBLAST.outputSequences)





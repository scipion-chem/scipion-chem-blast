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
from pyworkflow.protocol.params import PointerParam, BooleanParam, StringParam, EnumParam
from pyworkflow import BETA

from pwchem.utils.utilsFasta import fastFastaExport

from blast import Plugin, BLAST_DIC
from blast.constants import BLASTdbs

class ProtChemBLASTDatabase(EMProtocol):
    """Creates a BLAST database locally from a set of sequences or downloading from ncbi databases
    
    User IA Manual: BlastDatabase Protocol

The BlastDatabase protocol is designed to create local BLAST-compatible
databases from user-provided sequence files. This is a preparatory step
required for running local BLAST searches using the Blast protocol within
Scipion-Chem. It allows users to convert FASTA-formatted sequence collections
into searchable database formats that can be queried efficiently and repeatedly.

To begin, the user must provide a file containing one or more sequences in FASTA
format. These sequences can represent proteins or nucleotides, depending on the
intended use of the database. The protocol automatically detects the sequence
type, although the user can explicitly specify whether the database should be
created for protein or nucleotide queries. This selection determines how the
input is parsed and which search modes will be compatible with the resulting
database.

The user must also define a name for the database, which will be used to identify
it in downstream protocols. This name must be unique within the working project
and should not conflict with existing databases in the workspace. Optionally, a
title or description may be added to annotate the database and facilitate its
interpretation later on.

Once the protocol is executed, the input FASTA file is processed and indexed
using the `makeblastdb` utility. This step generates a collection of binary files
that store sequence data, indexing tables, and metadata required for fast
lookup. The output is stored in a designated folder and registered as a BLAST
database object within Scipion-Chem. This database can then be selected in the
Blast protocol to perform sequence similarity searches without requiring an
internet connection or access to remote servers.

In summary, the BlastDatabase protocol provides an efficient way to build
custom, local BLAST databases for personalized or large-scale sequence analysis.
It integrates seamlessly with the rest of the Scipion-Chem bioinformatics
toolkit and enables reproducible, offline workflows for sequence comparison and
homology-based modeling.
    """
    _label = 'Local BLAST database'
    _devStatus = BETA

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('fromNCBI', BooleanParam, default=False,
                      label='Create from a NCBI database: ')

        form.addParam('inputID', EnumParam, choices=BLASTdbs,
                      label='NCBI database name:', condition='fromNCBI',
                      help="NCBI database name to download. Options at 12/2021.\n"
                           "To get the current list use: {}".
                      format(os.path.join(Plugin.getVar(BLAST_DIC['home']), 'bin/update_blastdb.pl') + ' --showall'))

        group = form.addGroup('Input Sequences', condition='not fromNCBI')
        group.addParam('inputSequences', PointerParam, pointerClass='SetOfSequences',
                       label='Input sequences:', condition='not fromNCBI', allowsNull=True,
                       help="Set of sequences to create the database")
        group.addParam('dbType', EnumParam, default=0, condition='not fromNCBI',
                       choices=['Protein', 'Nucleotide'], display=EnumParam.DISPLAY_HLIST,
                       label='Type of database to query on: ')
        group.addParam('titleDB', StringParam,
                       label='New database name:', condition='not fromNCBI',
                       help="Name to designate the new database")


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.fromNCBI:
            self._insertFunctionStep('downloadDatabaseStep')
        else:
            self._insertFunctionStep('createDatabaseStep')

    def downloadDatabaseStep(self):
        dbName = self.getEnumText('inputID')
        args = ' --decompress {} -passive'.format(dbName)
        outDir = Plugin.getDatabasesDir()
        if not os.path.exists(outDir):
            os.mkdir(outDir)

        Plugin.updateDatabase(self, args, cwd=outDir)
        print('Database has been downloaded into {} directory'.format(outDir))

    def createDatabaseStep(self):
        inSeqs = self.inputSequences.get()
        inFasta = self._getPath('database.fasta')

        fastFastaExport(inSeqs, inFasta)
        dbClass = 'prot' if self.dbType.get() == 0 else 'nucl'
        outDir = Plugin.getDatabasesDir()

        args = ' -in {} -parse_seqids -title "{}" -dbtype {} -out {}'.\
          format(os.path.abspath(inFasta), self.titleDB.get(), dbClass, os.path.join(outDir, self.titleDB.get()))

        Plugin.runBLAST(self, 'makeblastdb', args, cwd=outDir)
        print('Database has been created as {} into {} directory'.format(self.titleDB.get(), outDir))


    def _validate(self):
        errors=[]
        return errors

    def _warnings(self):
        warns = []
        if not self.fromNCBI:
            if self.titleDB.get() in Plugin.getLocalDatabases():
                warns.append('There is already a database with that name in local.\n'
                             'If you continue, it will be overwritten.')
            elif self.titleDB.get() in BLASTdbs:
                warns.append('There is a NCBI database with that name.\n'
                             'If you continue, it may cause problems if you ever try to download it.')
        return warns

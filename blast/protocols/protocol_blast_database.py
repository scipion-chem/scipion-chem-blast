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
    
    AI Generated:

        ProtChemBLASTDatabase - User Manual

        Overview
        --------
        The ProtChemBLASTDatabase protocol allows users to create local BLAST-compatible
        databases from sequence collections. This is a preparatory step for running
        local BLAST searches using the ProtChemBLAST protocol within Scipion-Chem.
        Databases can be generated from user-provided FASTA sequences or downloaded
        directly from NCBI.

        Input Requirements
        ------------------
        1. **Create from NCBI**: Optionally, download a precompiled BLAST database
           from NCBI.
           - Specify the database name from the available NCBI BLASTdb list.
        2. **Local sequences**: Provide a `SetOfSequences` containing one or more
           sequences in FASTA format.
           - Specify whether sequences are **Protein** or **Nucleotide**.
           - Assign a **database name/title** for the new local database.

        Workflow
        --------
        1. **NCBI database download**:
           - Downloads and decompresses the selected database.
           - Stores files in the local BLAST database directory.
           - Ensures the database is ready for offline queries.

        2. **Local database creation**:
           - Exports the input sequences to a FASTA file.
           - Uses `makeblastdb` to create a BLAST database, generating
             necessary index and metadata files.
           - Supports both protein (`-dbtype prot`) and nucleotide (`-dbtype nucl`)
             databases.
           - The resulting database is registered in the local BLAST databases
             folder for use in ProtChemBLAST searches.

        Outputs
        -------
        - A local BLAST database ready for query by the ProtChemBLAST protocol.
        - Database files stored in the Scipion-Chem BLAST database directory.

        Validation & Warnings
        ---------------------
        - If a local database with the same name exists, a warning is issued
          to prevent accidental overwriting.
        - If the chosen name conflicts with an NCBI database, a warning
          is issued to avoid potential confusion in future downloads.
        - No additional numerical validation is required.

        Practical Recommendations
        -------------------------
        - Use NCBI downloads for standard reference databases.
        - Use local sequence input for custom or project-specific databases.
        - Ensure the database name is unique to prevent overwriting existing
          databases.
        - For nucleotide databases, ensure sequences are correctly formatted
          and free of ambiguous bases to avoid indexing errors.

        Final Perspective
        -----------------
        ProtChemBLASTDatabase enables reproducible and efficient preparation of
        BLAST databases, integrating seamlessly with Scipion-Chem workflows for
        sequence homology searches and annotation. It supports both offline
        custom databases and automated downloads from NCBI.
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

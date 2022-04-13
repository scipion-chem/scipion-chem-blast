# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

import pwem, os
from subprocess import check_call
from os.path import join, exists
from .constants import *

_version_ = '0.1'
_logo = "blast_logo.png"
_references = ['']

class Plugin(pwem.Plugin):
    _homeVar = BLAST_HOME
    _pathVars = [BLAST_HOME]
    _supportedVersions = [BLAST_DEFAULT_VERSION]

    _pluginHome = join(pwem.Config.EM_ROOT, BLAST + '-' + BLAST_DEFAULT_VERSION)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(BLAST_HOME, BLAST + '-' + BLAST_DEFAULT_VERSION)

    @classmethod
    def defineBinaries(cls, env):
        cls.addBLASTPackage(env)

    @classmethod
    def runBLAST(cls, protocol, program, args, cwd=None):
        """ Run BLAST program commands from a given protocol. """
        protocol.runJob(join(cls._pluginHome, 'bin', program), args, cwd=cwd)

    @classmethod
    def runEDirect(cls, protocol, fullCMD, cwd):
        """ Run EDirect program commands from a given protocol. """
        print('Edirect run: ', fullCMD)
        check_call(fullCMD, shell=True, env=protocol._getEnviron(), cwd=cwd)

    @classmethod
    def getEDirectProgram(cls, program):
        return join(cls._pluginHome, 'edirect', program)

    @classmethod
    def updateDatabase(cls, protocol, args, cwd=None):
        '''Upload a local BLAST database using the perl script'''
        if cwd is None:
            cwd = cls.getDatabasesDir()
        protocol.runJob('perl ' + join(cls._pluginHome, 'bin/update_blastdb.pl'), args, cwd=cwd)

    @classmethod  #  Test that
    def getEnviron(cls):
        pass


    @classmethod
    def addBLASTPackage(cls, env):
        installationCmd = 'wget %s -O %s && ' % (cls._getBLASTDownloadUrl(), cls._getBLASTTar())
        installationCmd += 'tar -xf %s --strip-components 1 && ' % cls._getBLASTTar()
        installationCmd += 'rm %s && ' % cls._getBLASTTar()
        installationCmd += 'mkdir %s && ' % join(pwem.Config.EM_ROOT, BLAST + '-' + BLAST_DEFAULT_VERSION, 'databases')

        #Edirect
        installationCmd += 'printf "N\\n" | sh -c "$(wget -q {} -O -)" && '.format(cls._getEDirectDownloadUrl())
        installationCmd += 'mv ~/edirect ./ && '

        # Creating validation file
        BLAST_INSTALLED = '%s_installed' % BLAST
        installationCmd += 'touch %s' % BLAST_INSTALLED  # Flag installation finished

        env.addPackage(BLAST,
                       version=BLAST_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=[(installationCmd, BLAST_INSTALLED)],
                       default=True)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def _getBLASTDownloadUrl(cls):
        return "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz"

    @classmethod
    def _getEDirectDownloadUrl(cls):
        return "ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh"

    @classmethod
    def _getBLASTTar(cls):
        pluginHome = join(pwem.Config.EM_ROOT, BLAST + '-' + BLAST_DEFAULT_VERSION)
        return pluginHome + '/' + BLAST + '-' + BLAST_DEFAULT_VERSION + '.tar.gz'

    @classmethod
    def getDatabasesDir(cls):
        return os.path.abspath(os.path.join(cls._pluginHome, 'databases'))

    @classmethod
    def getLocalDatabases(cls):
        databases = set([])
        for file in os.listdir(cls.getDatabasesDir()):
            databases.add(file.split('.')[0])
        databases = list(databases)
        databases.sort()
        return databases


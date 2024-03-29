# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

"""
Wizard to show default values for a blast configuration
"""

# Imports
from pwem.wizards.wizard import EmWizard
from pwchem.wizards import AddElementWizard
from .protocols import ProtChemBLAST, ProtChemNCBIDownload
from .constants import DEF_BLAST_PARAMS

class SetDefaultBLASTParameters(EmWizard):
    _targets = [(ProtChemBLAST, ['labelAdvice'])]

    def show(self, form):
        protocol = form.protocol
        form.setVar('evalue', 0.05)

        blatsProgram = protocol.getSelectedBLASTProgram()
        for attr in DEF_BLAST_PARAMS[blatsProgram]:
            form.setVar(attr, DEF_BLAST_PARAMS[blatsProgram][attr])


class AddNCBI_ID_Wizard(AddElementWizard):
    """Add ID or keyword in NCBI fetch protocol to the list"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
      inputParam, outputParam = self.getInputOutput(form)
      protocol = form.protocol

      inID = getattr(protocol, inputParam[0]).get()
      searchMode = getattr(protocol, inputParam[2]).get()

      if inID and inID.strip() != '':
          prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())

          if searchMode == 0:
              towrite = prevList + '{"ID": "%s"}\n' % (inID.strip())

          else:
              maxEntries = getattr(protocol, inputParam[1]).get()
              if maxEntries:
                  towrite = prevList + '{"ID": "%s", "maxEntries": "%s"}\n' % (inID.strip(), maxEntries)

          form.setVar(outputParam[0], towrite)


AddNCBI_ID_Wizard().addTarget(protocol=ProtChemNCBIDownload,
                              targets=['addEntry'],
                              inputs=['inputID', 'maxEntries', 'searchMode'],
                              outputs=['listIDs'])



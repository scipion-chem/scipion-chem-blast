================================
BLAST Scipion plugin
================================

Scipion framework plugin wrapping NCBI BLAST+ (2.12.0) for local sequence
similarity search, plus a set of protocols specific to the
B-Cell-Epitope-Prediction pipeline built on top of it:

- ``ProtChemBLAST``: generic local BLAST search.
- ``ProtChemBLASTDatabase``: local BLAST database construction.
- ``ProtChemNCBIDownload``: reference sequence retrieval from NCBI.
- ``ProtSelfToleranceFilter``: screens candidates against a local human
  proteome BLAST database (self-tolerance check), evaluated at the final
  candidate's own length rather than at its upstream parent region's.
- ``ProtBLASTPanelConservation``: scores each input ROI's conservation
  against a user-provided reference panel of strain/clade/variant
  sequences of the same pathogen under analysis.

BLAST+ itself is installed automatically (downloaded directly from the
NCBI FTP release). ``ProtSelfToleranceFilter`` additionally requires a
manually-built local human proteome BLAST database, pointed to via
``BLAST_HUMAN_DB_PATH`` in ``scipion.conf`` (building a full human-proteome
database is a heavyweight one-time step outside this protocol's own
scope).

===================
Install this plugin
===================

**Developer's version**

.. code-block::

            git clone https://github.com/Lvera-code/scipion-chem-blast.git
            cd scipion-chem-blast
            scipion3 installp -p . --devel

.. code-block::

            scipion3 tests blast.tests

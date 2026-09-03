=========
CHANGES
=========

0.1.0
=====
- Initial release: BLAST+ protocol suite for the B-Cell-Epitope-Prediction
  pipeline. ``ProtChemBLAST`` (generic BLAST search), ``ProtChemBLASTDatabase``
  (local database build), ``ProtChemNCBIDownload`` (reference sequence
  retrieval), ``ProtSelfToleranceFilter`` (human-proteome self-tolerance
  screening on the final candidate's own length, not its parent region's)
  and ``ProtBLASTPanelConservation`` (conservation scoring against a
  user-provided strain/clade/variant reference panel).

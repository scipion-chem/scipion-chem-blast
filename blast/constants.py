# **************************************************************************
# *
# * Authors: Daniel Del Hoyo Gomez
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


# Supported Versions
BLAST_DEFAULT_VERSION = '2.12.0'
BLAST, BLAST_HOME = 'blast', 'BLAST_HOME'

#BLAST PROTOCOL
MATCH, MATRIX, NOGAP = 0, 1, 2
dbProtChoices = ['Non-redundant (nr)', 'RefSeq Select (refseq_select_prot)', 'NCBI_Reference proteins (refseq_protein)',
                 'UniProtKB/Swiss-Prot (swissprot)', 'Patented protein (pataa)', 'Protein Data Bank (pdb)',
                 'Metagenomic (env_nr)', 'Transcriptome Shotgun Assembly (tsa_nr)',
                 'NCBI Mithocondrial Protein Reference Sequences (mito)',
                 'CDD database for delta-blast (cdd_delta)']

dbNucChoices = ['Nucleotide Collection (nt)', 'RefSeq Select RNA (refseq_select_rna)',
                'NCBI Reference RNA (refseq_rna)', 'RefSeq Representative genomes (refseq_representative_genomes)',
                'RefSeq Eukaryotic Representative Genomes (ref_euk_rep_genomes)',
                'RefSeq Prokaryote Representative Genomes (ref_prok_rep_genomes)',
                'Refseq viruses representative genomes (ref_viruses_rep_genomes)',
                'Refseq viroids representative genomes (ref_viroids_rep_genomes)',
                'RefSeq Genome Database (refseq_genomes)', 'Expressed sequence tags (est)',
                'Transcriptome Shotgun Assembly (tsa_nt)', 'Patent sequences (patnt)','PDB nucleotide database (pdbnt)',
                'Genomic survey sequences (gss)', 'Sequence tagged sites (dbsts)',
                'Environmental samples (env_nt)', 'NCBI Genomic Mithocondrial Reference Sequences (mito)',]

BLASTdbs = ['16S_ribosomal_RNA', '18S_fungal_sequences', '28S_fungal_sequences', 'Betacoronavirus', 'ITS_RefSeq_Fungi',
          'ITS_eukaryote_sequences', 'LSU_eukaryote_rRNA', 'LSU_prokaryote_rRNA', 'SSU_eukaryote_rRNA', 'env_nt',
          'env_nr', 'human_genome', 'landmark', 'mito', 'mouse_genome', 'nr', 'nt', 'pataa', 'patnt', 'pdbaa', 'pdbnt',
          'ref_euk_rep_genomes', 'ref_prok_rep_genomes', 'ref_viroids_rep_genomes', 'ref_viruses_rep_genomes',
          'refseq_select_rna', 'refseq_select_prot', 'refseq_protein', 'refseq_rna', 'swissprot', 'tsa_nr', 'tsa_nt',
          'taxdb']

matrixChoices = ['PAM30', 'PAM70', 'PAM250', 'BLOSUM80', 'BLOSUM62', 'BLOSUM45', 'BLOSUM50', 'BLOSUM90']

blastpProgramsHelp = 'BLASTP simply compares a protein query to a protein database.\n' \
                     'BLASTP-FAST is an accelerated version of BLASTP that is very fast and works best if the target ' \
                     'percent identity is 50% or more.\n' \
                     'PSI-BLAST allows the user to build a PSSM (position-specific scoring matrix) using the results ' \
                     'of the first BlastP run.\n' \
                     'DELTA-BLAST constructs a PSSM using the results of a Conserved Domain Database search and ' \
                     'searches a sequence database.'

blastnProgramsHelp = 'BLASTN is slow, but allows a word-size down to seven bases.\n' \
                     'Megablast is intended for comparing a query to closely related sequences and works best if the ' \
                     'target percent identity is 95% or more but is very fast.\n' \
                     'Discontiguous megablast uses an initial seed that ignores some bases (allowing mismatches) and ' \
                     'is intended for cross-species comparisons.'

#Configuration: match/mismatch: gap or matrix: gap
blastnMM = ['1/-2', '1/-3', '1/-4', '2/-3', '4/-5', '1/-1']
ALLOWED_BLAST_GAPS = {'1/-2': ['5/2', '2/2', '1/2', '0/2', '3/1', '2/1', '1/1'],
                      '1/-3': ['5/2', '2/2', '1/2', '0/2', '2/1', '1/1'],
                      '1/-4': ['5/2', '1/2', '0/2', '2/1', '1/1'],
                      '2/-3': ['4/4', '2/4', '0/4', '3/3', '6/2', '5/2', '4/2', '2/2'],
                      '4/-5': ['12/8', '6/5', '5/5', '4/5', '3/5'],
                      '1/-1': ['5/2', '3/2', '2/2', '1/2', '0/2', '4/1', '3/1', '2/1']}

ALLOWED_BLAST_GAPS.update({'PAM30': ['7/2', '6/2', '5/2', '10/1', '9/1', '8/1', '13/3', '15/3', '14/1', '14/2'],
                           'PAM70': ['8/2', '7/2', '6/2', '11/1', '10/1', '9/1', '12/3', '11/2'],
                           'PAM250': ['15/3', '14/3', '13/3', '12/3', '11/3', '17/2', '16/2', '15/2', '14/2', '21/1',
                                      '20/1', '19/1', '18/1', '17/1'],
                           'BLOSUM80': ['8/2', '7/2', '6/2', '11/1', '10/1', '9/1'],
                           'BLOSUM62': ['11/2', '10/2', '9/2', '8/2', '7/2', '6/2', '13/1', '12/1', '11/1', '10/1',
                                        '9/1'],
                           'BLOSUM45': ['13/3', '12/3', '11/3', '10/3', '15/2', '14/2', '13/2', '12/2', '19/1',
                                        '18/1', '17/1', '16/1'],
                           'BLOSUM50': ['13/3', '12/3', '11/3', '10/3', '9/3', '16/2', '15/2', '14/2', '13/2', '12/2',
                                        '19/1', '18/1', '17/1', '16/1', '15/1'],
                           'BLOSUM90': ['9/2', '8/2', '7/2', '6/2', '11/1', '10/1', '9/1']})

DEF_BLAST_PARAMS = {'blastn': {'word_size': '11', 'reward': '2', 'penalty': '-3', 'gapopen': '5', 'gapextend': '2'},
                    'megablast': {'word_size': '28', 'reward': '1', 'penalty': '-2',
                                  'gapopen': '', 'gapextend': ''},
                    'dc-megablast': {'word_size': '11', 'reward': '2', 'penalty': '-3',
                                     'gapopen': '5', 'gapextend': '2'},
                    'blastp': {'word_size': '6', 'matrix': 4, 'gapopen': '11', 'gapextend': '1'},
                    'blastp-fast': {'word_size': '6', 'matrix': 4, 'gapopen': '11', 'gapextend': '1'},
                    'psi-blast': {'word_size': '3', 'matrix': 4, 'gapopen': '11', 'gapextend': '1'},
                    'delta-blast': {'word_size': '3', 'matrix': 4, 'gapopen': '11', 'gapextend': '1'},
                    'blastx': {'word_size': '6', 'matrix': 4, 'gapopen': '11', 'gapextend': '1'},
                    'tblastn': {'word_size': '6', 'matrix': 4, 'gapopen': '11', 'gapextend': '1'},
                    'tblastx': {'word_size': '3', 'matrix': 4},
                    }







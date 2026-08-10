# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Enzo Sierra (enzogael57@gmail.com)
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
Shared batched-blastp-by-length-tier logic, used by both
ProtSelfToleranceFilter and ProtBLASTPanelConservation: ported from the
standalone B-Cell-Epitope-Prediction repo's src/engines/blast_engine.py
(_select_task/_select_evalue/_run_blastp_batch), which itself is reused
as-is (not duplicated) by that repo's own conservation_engine.py -- same
principle applied here, one copy shared by both protocols instead of two.

BLASTp's statistics depend heavily on query length: a strict default
E-value discards identical short-peptide hits as "not significant"; a lax
E-value on long queries produces irrelevant-homology noise. Peptides are
routed dynamically to a (task, E-value) tier by length, in batches per
tier (avoids re-invoking blastp once per peptide).
"""

import tempfile
from pathlib import Path
from typing import List, Tuple

import pandas as pd

OUTFMT6_COLUMNS: List[str] = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
]

SHORT_PEPTIDE_MAX_LEN = 30
MEDIUM_PEPTIDE_MAX_LEN = 100
EVALUE_SHORT = 50.0    # <= SHORT_PEPTIDE_MAX_LEN aa (lax: BLAST penalizes short queries)
EVALUE_MEDIUM = 0.1    # SHORT_PEPTIDE_MAX_LEN+1 .. MEDIUM_PEPTIDE_MAX_LEN aa
EVALUE_LONG = 0.05     # > MEDIUM_PEPTIDE_MAX_LEN aa (strict: avoids irrelevant-homology noise)


def select_task(sequence_length: int, short_max_len: int = SHORT_PEPTIDE_MAX_LEN) -> str:
    """'blastp-short' if sequence_length <= short_max_len, else 'blastp'."""
    return 'blastp-short' if sequence_length <= short_max_len else 'blastp'


def select_evalue(
    sequence_length: int,
    short_max_len: int = SHORT_PEPTIDE_MAX_LEN,
    medium_max_len: int = MEDIUM_PEPTIDE_MAX_LEN,
) -> float:
    """E-value tier by peptide length (see module docstring)."""
    if sequence_length <= short_max_len:
        return EVALUE_SHORT
    if sequence_length <= medium_max_len:
        return EVALUE_MEDIUM
    return EVALUE_LONG


def run_blastp_batch_via_runjob(protocol, plugin, records: List[Tuple[int, str]], task: str, db: str, evalue: float) -> pd.DataFrame:
    """Runs one length-tier batch of 'records' (idx, sequence) through blastp via the plugin's runBLAST/runJob.

    Args:
        protocol: The calling EMProtocol instance (needed by Plugin.runBLAST's runJob call).
        plugin: The blast Plugin class (has runBLAST classmethod).
        records: List of (original_index, sequence) -- all must share 'task'/'evalue'.
        task: 'blastp-short' or 'blastp'.
        db: BLAST database prefix (no extension).
        evalue: E-value for this tier.

    Returns:
        DataFrame in -outfmt 6 shape (OUTFMT6_COLUMNS), empty if 'records' is
        empty or blastp found no hits.
    """
    if not records:
        return pd.DataFrame(columns=OUTFMT6_COLUMNS)

    with tempfile.TemporaryDirectory(prefix=f'blastp_{task}_') as tmp:
        tmp_dir = Path(tmp)
        query_path = tmp_dir / 'query.fasta'
        out_path = tmp_dir / 'hits.tsv'

        with query_path.open('w', encoding='utf-8') as fh:
            for idx, seq in records:
                fh.write(f'>peptide_{idx}\n{seq}\n')

        args = (
            f'-task {task} -query {query_path} -db {db} -outfmt 6 '
            f'-evalue {evalue} -out {out_path}'
        )
        plugin.runBLAST(protocol, 'blastp', args)

        if not out_path.is_file() or out_path.stat().st_size == 0:
            return pd.DataFrame(columns=OUTFMT6_COLUMNS)
        return pd.read_csv(out_path, sep='\t', names=OUTFMT6_COLUMNS)


def run_batched_blastp(protocol, plugin, sequences: pd.Series, db: str) -> pd.DataFrame:
    """Runs 'sequences' (indexed by their original DataFrame index) through blastp, tiered by length.

    Returns:
        Concatenated -outfmt 6 hits across all tiers (OUTFMT6_COLUMNS), empty
        if 'sequences' is empty or no tier produced any hit.
    """
    if sequences.empty:
        return pd.DataFrame(columns=OUTFMT6_COLUMNS)

    tasks = sequences.str.len().apply(select_task)
    evalues = sequences.str.len().apply(select_evalue)

    frames = []
    tiers = pd.DataFrame({'task': tasks, 'evalue': evalues}).drop_duplicates()
    for task, evalue in tiers.itertuples(index=False):
        tier_mask = (tasks == task) & (evalues == evalue)
        records = list(zip(sequences.index[tier_mask], sequences[tier_mask]))
        frames.append(run_blastp_batch_via_runjob(protocol, plugin, records, task, db, evalue))

    non_empty = [f for f in frames if not f.empty]
    return pd.concat(non_empty, ignore_index=True) if non_empty else pd.DataFrame(columns=OUTFMT6_COLUMNS)


def max_identity_by_query(hits: pd.DataFrame, query_lengths: pd.Series, min_query_coverage: float) -> pd.Series:
    """Max %identity per query among hits covering >= min_query_coverage of the query's own length.

    Without the coverage filter, a tiny 5-6aa fragment 100% identical by
    pure chance against a large database would reject peptides that don't
    actually resemble anything real (ported from
    blast_engine._max_identity_by_query).
    """
    if hits.empty:
        return pd.Series(dtype=float)

    hit_query_idx = hits['qseqid'].str.replace('peptide_', '', regex=False).astype(int)
    hit_query_length = hit_query_idx.map(query_lengths)
    coverage = hits['length'] / hit_query_length
    covered = hits[coverage >= min_query_coverage]

    if covered.empty:
        return pd.Series(dtype=float)
    return covered.groupby('qseqid')['pident'].max()


def panel_breadth_by_query(hits: pd.DataFrame, query_lengths: pd.Series, identity_threshold: float, min_query_coverage: float) -> pd.Series:
    """Count of DISTINCT panel sequences each query matches above the thresholds (breadth, not best-hit).

    Ported from conservation_engine._panel_breadth_by_query: measures how
    many different panel members (sseqid) a candidate resembles, not just
    the single best match -- see ProtBLASTPanelConservation for why breadth
    (not best-hit) is the right metric for cross-strain conservation.
    """
    if hits.empty:
        return pd.Series(dtype=int)

    hit_query_idx = hits['qseqid'].str.replace('peptide_', '', regex=False).astype(int)
    hit_query_length = hit_query_idx.map(query_lengths)
    coverage = hits['length'] / hit_query_length
    qualifying = hits[(coverage >= min_query_coverage) & (hits['pident'] >= identity_threshold)]

    if qualifying.empty:
        return pd.Series(dtype=int)
    return qualifying.groupby('qseqid')['sseqid'].nunique()

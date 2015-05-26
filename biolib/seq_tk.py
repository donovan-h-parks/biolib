###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'


import string


"""Sequence manipulation and statistics."""

_complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')


def count_nt(seq):
    """Count occurences of each nucleotide in a sequence.

    Only the bases A, C, G, and T(U) are counted. Ambiguous
    bases are ignored.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    list
        Number of A, C, G, and T(U) in sequence.
    """

    s = seq.upper()
    a = s.count('A')
    c = s.count('C')
    g = s.count('G')
    t = s.count('T') + s.count('U')

    return a, c, g, t


def gc(seq):
    """Calculate GC content of a sequence.

    GC is calculated as (G+C)/(A+C+G+T), where
    each of these terms represents the number
    of nucleotides within the sequence. Ambiguous
    and degenerate bases are ignored. Uracil (U)
    is treated as a thymine (T).

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    float
        GC content of sequence.
    """

    a, c, g, t = count_nt(seq)
    return float(g + c) / (a + c + g + t)


def ambiguous_nucleotides(seq):
    """Count ambiguous or degenerate nucleotides in a sequence.

    Any base that is not a A, C, G, or T/U is considered
    to be ambiguous or degenerate.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    int
        Number of ambiguous and degenerate bases.
    """

    a, c, g, t = count_nt(seq)
    return len(seq) - (a + c + g + t)


def rev_comp(seq):
    """Rverse complement a sequence."""
    return seq.translate(_complements)[::-1]


def N50(seqs):
    """Calculate N50 for a set of sequences.

     N50 is defined as the length of the longest
     sequence, L, for which 50% of the total bases
     are present in sequences of length >= L.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    int
        N50 for the set of sequences.
    """

    seq_lens = [len(x) for x in seqs.values()]
    threshold = sum(seq_lens) / 2.0

    seq_lens.sort(reverse=True)

    current_sum = 0
    for seq_len in seq_lens:
        current_sum += seq_len
        if current_sum >= threshold:
            N50 = seq_len
            break

    return N50


def mean_length(seqs):
    """Calculate mean length of sequences.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    int
        Mean length of sequences.
    """

    total_len = 0
    for _seq_id, seq in seqs:
        total_len += len(seq)

    return float(total_len) / len(seqs)


def max_length(seqs):
    """Identify longest sequence.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence ids.

    Returns
    -------
    int
        Length of longest sequence.
    """

    longest_seq = 0
    for _seq_id, seq in seqs:
        longest_seq = max(longest_seq, len(seq))

    return longest_seq


def fragment(seq, window_size, step_size):
    """Fragment sequence into fixed sized windows.

    The last fragment may not be shorter than
    the window size, but will only be generated
    if it is at least half the window size.

    Parameters
    ----------
    seq : str
        Sequence to fragment.
    window_size : int
        Size of each fragment.
    step_size : int
        Number of bases to move after each window.

    Returns
    -------
    list
        Fragments from sequences.
    """

    fragments = []
    start = 0
    for i in xrange(0, len(seq), step_size):
        end = i + window_size
        if end < len(seq):
            fragments.append(seq[start:end])
            start = end

    # get last fragment if it is at least half
    # the specified window size
    if len(seq) - start >= 0.5 * window_size:
        fragments.append(seq[start:])

    return fragments


def trim_seqs(seqs, min_per_taxa, min_per_bp):
    """Trim multiple sequence alignment.

    Parameters
    ----------
    seqs : d[seq_id] -> sequence
        Aligned sequences.
    min_per_taxa : float
        Minimum percentage of taxa required to retain a column.
    min_per_bp : float
        Minimum percentage of base pairs required to keep trimmed sequence.

    Returns
    -------
    dict : d[seq_id] -> sequence
        Dictionary of trimmed sequences.
    dict : d[seq_id] -> sequence
        Dictionary of pruned sequences.
    """

    alignment_length = len(seqs.values()[0])

    # count number of taxa represented in each column
    column_count = [0] * alignment_length
    for seq in seqs.values():
        for i, ch in enumerate(seq):
            if ch != '.' and ch != '-':
                column_count[i] += 1

    mask = [False] * alignment_length
    for i, count in enumerate(column_count):
        if count >= min_per_taxa * len(seqs):
            mask[i] = True

    # trim leading and trailing columns
    output_seqs = {}
    pruned_seqs = {}
    for seq_id, seq in seqs.iteritems():
        masked_seq = ''.join([seq[i] for i in xrange(0, len(mask)) if mask[i]])

        valid_bases = len(masked_seq) - masked_seq.count('.') - masked_seq.count('-')
        if valid_bases < len(masked_seq) * min_per_bp:
            pruned_seqs[seq_id] = masked_seq
            continue

        output_seqs[seq_id] = masked_seq

    return output_seqs, pruned_seqs

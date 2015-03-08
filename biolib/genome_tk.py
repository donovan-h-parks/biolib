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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

from collections import defaultdict

from biolib.common import remove_extension
from biolib.seq_io import read_fasta_seq

"""
Functions for verify, exploring, modifying and
calculating statistics on one or more genomes.
"""


def unique(genome_files):
    """Check if sequences are assigned to multiple bins.

    Parameters
    ----------
    genome_files : iterable
        Path to genome fasta files.

    Returns
    -------
    dict : d[genome_id][genome_id] -> [shared sequences]
        List of any sequences within a genome observed multiple times.
    """

    # read sequence IDs from all genomes,
    # while checking for duplicate sequences within a genomes
    duplicates = defaultdict(lambda: defaultdict(list))

    genome_seqs = {}
    for f in genome_files:
        genome_id = remove_extension(f)

        seq_ids = set()
        for seq_id, _seq in read_fasta_seq(f):
            if seq_id in seq_ids:
                duplicates[genome_id][genome_id].append(seq_id)

            seq_ids.add(seq_id)

        genome_seqs[genome_id] = seq_ids

    # check for sequences assigned to multiple bins
    genome_ids = genome_seqs.keys()
    for i in xrange(0, len(genome_ids)):
        seq_idsI = genome_seqs[genome_ids[i]]

        for j in xrange(i + 1, len(genome_ids)):
            seq_idsJ = genome_seqs[genome_ids[j]]

            seq_intersection = seq_idsI.intersection(seq_idsJ)

            if len(seq_intersection) > 0:
                duplicates[genome_ids[i]][genome_ids[j]] = seq_intersection
                duplicates[genome_ids[j]][genome_ids[i]] = seq_intersection

    return duplicates

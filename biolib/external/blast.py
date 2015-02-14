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

import os
import sys
import logging
import subprocess


class BlastRunner():
    """Wrapper for running blast."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        self._check_for_blast()

    def blastp(self, query_seqs, prot_db, evalue, cpus, output_file):
        """Apply blastp to query file.

        Finds homologs to query sequences using blastp homology search
        against a protein database.

        Parameters
        ----------
        query_seqs : str
            File containing query sequences.
        prot_db : str
            File containing blastp formatted database.
        evalue : float
            E-value threshold used to identify homologs.
        cpus : int
            Number of cpus to use during homology search.
        output_file : str
            Output file containing blastp results.
        """

        cmd = "blastp -num_threads %d" % cpus
        cmd += " -query %s -db %s -out %s -evalue %g" % (query_seqs, prot_db, output_file, evalue)
        cmd += " -outfmt '6 qseqid qlen sseqid slen length pident evalue bitscore'"
        os.system(cmd)

    def _check_for_blast(self):
        """Check to see if blastp is on the system before we try to run it."""

        # Assume that a successful blast -help returns 0 and anything
        # else returns non-zero
        try:
            subprocess.call(['blastp', '-help'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  [Error] Make sure blastp is on your system path.")
            sys.exit()


class BlastParser():
    """Parses output files produced with BlastRunner."""

    def __init__(self):
        """Initialization."""
        pass

    def identify_homologs(self, blast_table, evalue_threshold, per_identity_threshold, per_aln_len_threshold):
        """Identify homologs among  blast hits.

        Identifies hits satisfying the criteria required for a
        gene to be considered a homolog.

        Parameters
        ----------
        blast_table : str
            File containing blast hits in the custom tabular format produced by BlastRunner.
        evalue_threshold : float
            E-value threshold used to define homologous gene.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.

        Returns
        -------
        set
            Identifiers for homologous genes.
        """

        homologs = set()
        for line in open(blast_table):
            line_split = line.split('\t')

            _query_seq_id = line_split[0]
            query_len = int(line_split[1])

            sub_seq_id = line_split[2]
            _subject_len = int(line_split[3])

            aln_len = int(line_split[4])
            per_ident = float(line_split[5])
            evalue = float(line_split[6])
            _bitscore = float(line_split[7])

            if evalue <= evalue_threshold and per_ident >= per_identity_threshold:
                per_aln_len = aln_len * 100.0 / query_len

                if per_aln_len >= per_aln_len_threshold:
                    homologs.add(sub_seq_id)

        return homologs

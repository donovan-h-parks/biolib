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
import logging

from biolib.external.execute import check_on_path

"""
To do:
 - this should be extended to handle blastn and some of the variants
 - this class and the diamond class should mirror each other
   to the extent possible
"""


class Blast():
    """Wrapper for running blast."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger()

        check_on_path('blastp')

        self.cpus = cpus

        self.output_fmt = {'standard': '6',
                            'custom': '6 qseqid qlen sseqid slen length pident evalue bitscore'}

    def blastp(self, query_seqs, prot_db, evalue, output_fmt, output_file):
        """Apply blastp to query file.

        Finds homologs to query sequences using blastp homology search
        against a protein database. Hit can be reported using  either
        the 'standard' table 6 format or the following 'custom' format:
            qseqid qlen sseqid slen length pident evalue bitscore


        Parameters
        ----------
        query_seqs : str
            File containing query sequences.
        prot_db : str
            File containing blastp formatted database.
        evalue : float
            E-value threshold used to identify homologs.
        output_fmt : str
            Specified output format of blast table: standard or custom.
        output_file : str
            Output file containing blastp results.
        """

        cmd = "blastp -num_threads %d" % self.cpus
        cmd += " -query %s -db %s -out %s -evalue %g" % (query_seqs, prot_db, output_file, evalue)
        cmd += " -outfmt %s" % self.output_fmt[output_fmt]
        os.system(cmd)

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

import os
import tempfile
import logging

from biolib.external.execute import check_on_path


class Diamond(object):
    """Wrapper for running diamond."""

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger()

        check_on_path('diamond')

        self.cpus = cpus

    def make_database(self, prot_file, db_file):
        """Make diamond database.

        Parameters
        ----------
        prot_file : str
            Fasta file with protein sequences.
        db_file : str
            Desired name of Diamond database.
        """

        cmd = 'diamond makedb -p %d --in %s -d %s' % (self.cpus, prot_file, db_file)
        os.system(cmd)

    def blastp(self, prot_file, db_file, evalue, per_identity, max_target_seqs, diamond_daa_file):
        """Apply diamond blastp to a set of protein sequences.

        Parameters
        ----------
        prot_file : str
            Fasta file with protein sequences.
        db_file : str
            Diamond database of protein sequences.
        evalue : float
            E-value threshold used by blast.
        per_identity : float
            Percent identity threshold used by blast [0, 100].
        max_target_seqs : int
            Maximum number of hits to report per sequence.
        diamond_daa_file : str
            Desired name of Diamond data file.
        """

        if db_file.endswith('.dmnd'):
            db_file = db_file[0:db_file.rfind('.dmnd')]

        cmd = "diamond blastp --seg no -p %d -t %s -q %s -d %s -e %g --id %f -k %d -a %s" % (self.cpus,
                                                                                tempfile.gettempdir(),
                                                                                prot_file,
                                                                                db_file,
                                                                                evalue,
                                                                                per_identity,
                                                                                max_target_seqs,
                                                                                diamond_daa_file)

        os.system(cmd)

    def blastx(self, nt_file, db_file, evalue, per_identity, max_target_seqs, diamond_daa_file):
        """Apply diamond blastx to a set of nucleotide sequences.

        Parameters
        ----------
        nt_file : str
            Fasta file with nucleotide sequences.
        db_file : str
            Diamond database of protein sequences.
        evalue : float
            E-value threshold used by blast.
        per_identity : float
            Percent identity threshold used by blast [0, 100].
        max_target_seqs : int
            Maximum number of hits to report per sequence.
        diamond_daa_file : str
            Desired name of Diamond data file.
        """

        if db_file.endswith('.dmnd'):
            db_file = db_file[0:db_file.rfind('.dmnd')]

        cmd = 'diamond blastx --seg no -p %d -t %s -q %s -d %s -e %f --id %f -k %d -a %s' % (self.cpus,
                                                                                        tempfile.gettempdir(),
                                                                                        nt_file,
                                                                                        db_file,
                                                                                        evalue,
                                                                                        per_identity,
                                                                                        max_target_seqs,
                                                                                        diamond_daa_file)

        os.system(cmd)

    def view(self, diamond_daa_file, output_table, compress=False):
        """Generate flat file from diamond DAA file.

        Parameters
        ----------
        diamond_daa_file : str
            Diamond DAA file.
        output_table : str
            Diamond database of protein sequeces.
        compress : boolean
            Flag indicating if output table should be compressed.
        """

        cmd = 'diamond view -p %d -a %s -o %s --compress %d' % (self.cpus,
                                                                    diamond_daa_file,
                                                                    output_table,
                                                                    compress)
        os.system(cmd)

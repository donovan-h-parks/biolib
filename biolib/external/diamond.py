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
import sys
import subprocess
import logging

"""
To do:
 - this class and the blast class should mirror each other
   to the extent possible
"""


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

        self._check_for_diamond()

        self.cpus = cpus

    def _check_for_diamond(self):
        """Check to see if BLAST is on the system path."""
        try:
            subprocess.call(['diamond', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("[Error] Make sure diamond is on your system path.")
            sys.exit(-1)

    def blastx(self, nt_file, db_file, evalue, per_identity, max_target_seqs, output_file):
        """Apply diamond blastx to a set of nucleotide sequences.

        Parameters
        ----------
        nt_file : str
            Fasta file with nucleotide sequences.
        db_file : str
            Diamond database of protein sequeces.
        evalue : float
            E-value threshold used by blast.
        per_identity : float
            Percent identity threshold used by blast.
        max_target_seqs : int
            Maximum number of hits to report per sequence.
        output_file : str
            File to store hits identified by diamond.
        """

        if db_file.endswith('.dmnd'):
            db_file = db_file[0:db_file.rfind('.dmnd')]

        os.system('diamond blastx --compress 0 -p %d -q %s -d %s -e %f --id %f -k %d -o %s' % (self.cpus,
                                                                            nt_file,
                                                                            db_file,
                                                                            evalue,
                                                                            per_identity,
                                                                            max_target_seqs,
                                                                            output_file))

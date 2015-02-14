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


class MuscleError(BaseException):
    pass


class MuscleRunner():
    """Wrapper for running muscle."""

    def __init__(self):
        """Initialization."""
        self.logger = logging.getLogger()

        self._check_for_muscle()

    def run(self, seqs, output_file, log_file):
        """Apply muscle to query sequences.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequence to alignment.
        output_file: str
            Output file containing multiple sequence alignment.
        log_file: str
            Output file containing information about running of muscle.
        """

        cmd = 'muscle -quiet -in %s -out %s -log %s' % (seqs, output_file, log_file)
        os.system(cmd)

    def _check_for_muscle(self):
        """Check to see if muscle is on the system before we try to run it."""

        # Assume that a successful blast -help returns 0 and anything
        # else returns non-zero
        try:
            subprocess.call(['muscle', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  [Error] Make sure muscle is on your system path.")
            sys.exit()

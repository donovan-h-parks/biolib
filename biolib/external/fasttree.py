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

import os
import sys
import logging

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"


class FastTreeRunner():
    """Wrapper for running FastTree."""

    def __init__(self, multithreaded=True):
        """Initialization."""
        self.logger = logging.getLogger()

        self._check_for_fasttree()

        self.multithreaded = multithreaded

    def run(self, msa_file, model_str, output_tree, output_tree_log, log_file=None):
        """Infer tree using FastTree.

        Parameters
        ----------
        msa_file : str
            Fasta file containing multiple sequence alignment.
        model_str : str
            Specified either the 'wag' or 'jtt' model.
        output_tree: str
            Output file containing inferred tree.
        output_tree_log: str
            Output file containing information about inferred tree.
        output_log: str
            Output file containing information about running of FastTree.
        """

        if model_str.upper() == 'JTT':
            model_str = ''
        elif model_str.upper() == 'WAG':
            model_str = '-wag'

        if not log_file:
            log_file = '/dev/null'

        cmd = '-quiet -nosupport -gamma %s -log %s %s > %s 2> %s' % (model_str, output_tree_log, msa_file, output_tree, log_file)
        if self.multithreaded:
            cmd = 'FastTreeMP ' + cmd
            os.system(cmd)
        else:
            cmd = 'FastTree ' + cmd
            os.system(cmd)

    def _check_for_fasttree(self):
        """Check to see if FastTree is on the system path."""

        try:
            exit_status = os.system('FastTree 2> /dev/null')
        except:
            print "Unexpected error!", sys.exc_info()[0]
            raise

        if exit_status != 0:
            print "[Error] FastTree is not on the system path."
            sys.exit()

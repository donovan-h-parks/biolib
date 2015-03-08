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
import logging

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

from biolib.external.execute import check_on_path


class FastTree():
    """Wrapper for running FastTree."""

    def __init__(self, multithreaded=True):
        """Initialization."""
        self.logger = logging.getLogger()

        self.multithreaded = multithreaded

        if self.multithreaded:
            check_on_path('FastTreeMP')
        else:
            check_on_path('FastTree')

    def run(self, msa_file, seq_type, model_str, output_tree, output_tree_log, log_file=None):
        """Infer tree using FastTree.

        Parameters
        ----------
        msa_file : str
            Fasta file containing multiple sequence alignment.
        seq_type : str
            Specifies multiple sequences alignment is of 'nt' or 'prot'.
        model_str : str
            Specified either the 'wag' or 'jtt' model.
        output_tree: str
            Output file containing inferred tree.
        output_tree_log: str
            Output file containing information about inferred tree.
        output_log: str
            Output file containing information about running of FastTree.
        """

        if seq_type.upper() == 'PROT':
            seq_type_str = ''
        elif seq_type.upper() == 'NT':
            seq_type_str = '-nt'

        if model_str.upper() == 'JTT':
            model_str = ''
        elif model_str.upper() == 'WAG':
            model_str = '-wag'

        if not log_file:
            log_file = '/dev/null'

        cmd = '-quiet -nosupport -gamma %s %s -log %s %s > %s 2> %s' % (seq_type_str,
                                                                        model_str,
                                                                        output_tree_log,
                                                                        msa_file,
                                                                        output_tree,
                                                                        log_file)
        if self.multithreaded:
            cmd = 'FastTreeMP ' + cmd
            os.system(cmd)
        else:
            cmd = 'FastTree ' + cmd
            os.system(cmd)

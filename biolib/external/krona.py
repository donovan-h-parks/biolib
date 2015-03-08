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

""" This is just a stub right now. """


class Krona():
    """Wrapper for creating Krona plots with KronaTools."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()

        check_on_path('ktImportText')

    def create(self):
        """Create Krona plot from ...

        """
        pass
        #cmd = ["ktImportText",'-o',outputName]
        #for i, tmp in enumerate(tempfile_paths):
        #    cmd.append(','.join([tmp,otuTables[i].sample_name]))

        # run the actual krona
        #subprocess.check_call(' '.join(cmd), shell=True)

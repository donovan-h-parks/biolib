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
import ntpath
from StringIO import StringIO

from biolib.common import make_sure_path_exists

def logger_setup(log_dir, log_file, program_name, version, silent):
    """Set logging for application.

    Parameters
    ----------
    log_dir : str
        Output directory for log file.
    log_file : str
        Desired name of log file.
    program_name : str
        Name of program.
    version : str
        Program version number.
    silent : boolean
        Flag indicating if output to stdout should be suppressed.
    """

    # setup general properties of logger
    logger = logging.getLogger('')
    logger.setLevel(logging.DEBUG)
    log_format = logging.Formatter(fmt="[%(asctime)s] %(levelname)s: %(message)s",
                                   datefmt="%Y-%m-%d %H:%M:%S")

    # setup logging to console
    stream_logger = logging.StreamHandler(sys.stdout)
    stream_logger.setFormatter(log_format)
    stream_logger.setLevel(logging.DEBUG)
    logger.addHandler(stream_logger)
    if silent:
        stream_logger.setLevel(logging.ERROR)
        sys.stdout = StringIO()

    if log_dir:
        make_sure_path_exists(log_dir)
        file_logger = logging.FileHandler(os.path.join(log_dir, log_file), 'a')
        file_logger.setFormatter(log_format)
        logger.addHandler(file_logger)

    logger.info('%s v%s' % (program_name, version))
    logger.info(ntpath.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:]))
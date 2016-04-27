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
import shutil
import tempfile
import ntpath

import biolib.seq_io as seq_io
from biolib.parallel import Parallel
from biolib.bootstrap import bootstrap_alignment, bootstrap_support
from biolib.external.execute import check_on_path


class FastTree():
    """Wrapper for running FastTree."""

    def __init__(self, multithreaded=True):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.multithreaded = multithreaded

        if self.multithreaded:
            check_on_path('FastTreeMP')
        else:
            check_on_path('FastTree')

    def bootstrap(self, input_tree, msa_file, seq_type, model_str, num_replicates, output_tree, cpus):
        """Perform non-parametric bootstrapping.

        Parameters
        ----------
        input_tree : str
            File containing newick tree to decorate with bootstraps.
        msa_file : str
            Fasta file containing multiple sequence alignment.
        seq_type : str
            Specifies multiple sequences alignment is of 'nt' or 'prot'.
        model_str : str
            Specified either the 'wag' or 'jtt' model.
        num_replicates : int
            Number of replicates to perform.
        output_tree: str
            Output file containing tree with bootstrap values.
        cpus : int
            Number of cpus to use.
        """

        assert(seq_type in ['nt', 'prot'])
        assert(model_str in ['wag', 'jtt'])

        self.replicate_dir = tempfile.mkdtemp()
        self.seq_type = seq_type
        self.model = model_str
        self.msa = seq_io.read(msa_file)

        # calculate replicates
        parallel = Parallel(cpus)
        parallel.run(self._bootstrap, None, xrange(num_replicates), None)

        # calculate support values
        rep_tree_files = []
        for rep_index in xrange(num_replicates):
            rep_tree_files.append(os.path.join(self.replicate_dir, 'bootstrap.tree.' + str(rep_index) + '.tre'))

        bootstrap_support(input_tree, rep_tree_files, output_tree)

        shutil.rmtree(self.replicate_dir)

    def _bootstrap(self, replicated_num):
        """Infer tree from bootstrapped multiple sequence alignment.

        Parameters
        ----------
        replicated_num : int
          Unique replicate number.
        cpus : int
            Number of cpus to use.
        """

        output_msa = os.path.join(self.replicate_dir, 'bootstrap.msa.' + str(replicated_num) + '.fna')
        bootstrap_alignment(self.msa, output_msa)

        output_tree = os.path.join(self.replicate_dir, 'bootstrap.tree.' + str(replicated_num) + '.tre')
        fast_tree_output = os.path.join(self.replicate_dir, 'bootstrap.fasttree.' + str(replicated_num) + '.out')
        self.run(output_msa, self.seq_type, self.model, output_tree, fast_tree_output)

        return True

    def parallel_run(self, msa_files, seq_type, model_str, output_dir, cpus):
        """Infer tree using FastTree in parallel.

        Parameters
        ----------
        msa_files : str
            Fasta files containing multiple sequence alignments.
        seq_type : str
            Specifies multiple sequences alignment is of 'nt' or 'prot'.
        model_str : str
            Specified either the 'wag' or 'jtt' model.
        output_dir: str
            Prefix for all output files.
        """

        assert(seq_type in ['nt', 'prot'])
        assert(model_str in ['wag', 'jtt'])

        self.output_dir = output_dir
        self.seq_type = seq_type
        self.model = model_str

        parallel = Parallel(cpus)
        parallel.run(self._parallel_infer_tree, None, msa_files, None)

    def _parallel_infer_tree(self, msa_file):
        """Infer tree using FastTree in parallel.

        Parameters
        ----------
        msa_file : str
            Fasta files containing multiple sequence alignments.
        """

        file_prefix = ntpath.basename(msa_file)
        if '.' in file_prefix:
            file_prefix = file_prefix[0:file_prefix.find('.')]

        output_tree = os.path.join(self.output_dir, file_prefix + '.tree')
        fast_tree_output = os.path.join(self.output_dir, file_prefix + '.log')
        self.run(msa_file, self.seq_type, self.model, output_tree, fast_tree_output)

    def run(self, msa_file, seq_type, model_str, output_tree, output_tree_log, log_file=None):
        """Infer tree using FastTree.

        All trees are inferred using the GAMMA distribution to model
        rate heterogeneity. Nucleotide trees are inferred under the
        GTR model.

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

        assert(seq_type in ['nt', 'prot'])
        assert(model_str in ['wag', 'jtt'])

        if seq_type == 'prot':
            seq_type_str = ''
            if model_str.upper() == 'JTT':
                model_str = ''
            elif model_str.upper() == 'WAG':
                model_str = '-wag'
        elif seq_type == 'nt':
            seq_type_str = '-nt'
            model_str = '-gtr'

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

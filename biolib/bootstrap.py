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

import logging
import random

import dendropy


class Bootstrap(object):
    """Perform non-parametric bootstrapping on multiple sequence alignment."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()

    def support_values(self, input_tree, replicate_trees, output_tree):
        """ Calculate support for tree with replicates covering the same taxon set.

        Parameters
        ----------
        input_tree : str
          Tree inferred from complete data.
        replicate_trees : iterable
          Files containing replicate trees.
        output_tree: str
          Name of output tree with support values.
        """

        tree = dendropy.Tree.get_from_path(input_tree, schema='newick', rooting="force-unrooted", preserve_underscores=True)

        rep_trees = []
        for rep_tree_file in replicate_trees:
            rep_trees.append(dendropy.Tree.get_from_path(rep_tree_file, schema='newick', rooting="force-unrooted", preserve_underscores=True))

        rep_tree_list = dendropy.TreeList(rep_trees)

        for node in tree.internal_nodes():
            taxa_labels = [x.taxon.label for x in node.leaf_nodes()]
            bootstrap = int(rep_tree_list.frequency_of_bipartition(labels=taxa_labels) * 100)

            if node.label:
                node.label = str(bootstrap) + ':' + node.label
            else:
                node.label = str(bootstrap)

        tree.write_to_path(output_tree, schema='newick', suppress_rooting=True, unquoted_underscores=True)

    def bootstrap(self, msa, output_file):
        """Bootstrap multiple sequence alignment.

        Parameters
        ----------
        msa : d[seq_id] -> seq
          Full multiple sequence alignment.
        output_file : str
          File to write bootstrapped alignment.
        """
        alignment_len = len(msa[msa.keys()[0]])
        cols = [random.randint(0, alignment_len - 1) for _ in xrange(alignment_len)]

        fout = open(output_file, 'w')
        for seq_id, seq in msa.iteritems():
            fout.write('>' + seq_id + '\n')
            for col in cols:
                fout.write(seq[col])
            fout.write('\n')
        fout.close()

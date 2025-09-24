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

"""Perform non-parametric bootstrapping on multiple sequence alignment."""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'


import random
import logging
from multiprocessing import Pool, cpu_count
from collections import Counter

import toytree
from tqdm import tqdm


def get_bipartitions_from_tree(tree_path):
    """Extract bipartitions from a tree using toytree."""

    try:
        # Load tree with toytree
        tree = toytree.tree(tree_path)  
        tree = tree.unroot()
        
        # Get all tip names
        tip_names = tree.get_tip_labels()
        tip_set = set(tip_names)
        
        bipartitions = []
        
        # Get internal nodes (exclude tips and root)
        internal_nodes = [node for node in tree.treenode.traverse() 
                         if not node.is_leaf() and not node.is_root()]
        
        for node in internal_nodes:
            # Get all tip names below this node
            descendants = set()
            for leaf in node.get_leaves():
                if hasattr(leaf, 'name') and leaf.name:
                    descendants.add(leaf.name)
                elif hasattr(leaf, 'idx'):
                    # Use tip label by index if name not available
                    descendants.add(tip_names[leaf.idx])
            
            # Create bipartition - use smaller set for consistent representation
            if len(descendants) > 0 and len(descendants) < len(tip_set):
                if len(descendants) <= len(tip_set) / 2:
                    bipartition = frozenset(descendants)
                else:
                    bipartition = frozenset(tip_set - descendants)
                
                bipartitions.append(bipartition)
        
        return bipartitions, tip_names
        
    except Exception as e:
        logging.getLogger().error(f"Error processing tree {tree_path}: {e}")
        return [], []


def process_tree(args):
    """Worker: process a single replicate tree file and return split counts."""

    tree_path, reference_tips = args
    
    bipartitions, _tip_names = get_bipartitions_from_tree(tree_path)
    
    # Filter bipartitions to only include those consistent with reference
    valid_bipartitions = []
    reference_set = set(reference_tips)
    
    for bp in bipartitions:
        # Only keep bipartitions where all tips are in reference set
        if bp.issubset(reference_set):
            valid_bipartitions.append(bp)
    
    return Counter(valid_bipartitions)


def apply_support(tree, split_counts, total_trees):
    """Apply bootstrap support values to toytree object."""
    
    tip_names = tree.get_tip_labels()
    tip_set = set(tip_names)
    
    # Get internal nodes
    internal_nodes = [node for node in tree.treenode.traverse() 
                     if not node.is_leaf() and not node.is_root()]

    for node in internal_nodes:
        # Get descendants for this node
        descendants = set()
        for leaf in node.get_leaves():
            if hasattr(leaf, 'name') and leaf.name:
                descendants.add(leaf.name)
            elif hasattr(leaf, 'idx'):
                descendants.add(tip_names[leaf.idx])
        
        if len(descendants) > 0 and len(descendants) < len(tip_set):
            # Create bipartition
            if len(descendants) <= len(tip_set) / 2:
                bipartition = frozenset(descendants)
            else:
                bipartition = frozenset(tip_set - descendants)
            
            # Get support value
            support_count = split_counts.get(bipartition, 0)
            support_percentage = (support_count / total_trees) * 100
            
            # Set support value on node
            node.support = support_percentage


def bootstrap_support(input_tree, replicate_trees, output_tree, cpus=None):
    """Determine bootstrap support for tree using toytree."""
    
    logger = logging.getLogger()
    
    # Load original tree to get reference tip names
    logger.info('Reading original tree.')
    orig_tree = toytree.tree(input_tree)

    # Get children in left and right branches before unrooting
    left_root_tips = None
    right_root_tips = None
    if len(orig_tree.treenode.children) == 2:  # Only for rooted trees
        left_child = orig_tree.treenode.children[0]
        right_child = orig_tree.treenode.children[1]
        
        # Get tip names in left and right subtrees
        left_root_tips = set([leaf.name for leaf in left_child.get_leaves()])
        right_root_tips = set([leaf.name for leaf in right_child.get_leaves()])

    orig_tree = orig_tree.unroot()
    reference_tips = orig_tree.get_tip_labels()
    
    if cpus is None:
        cpus = cpu_count()
    
    # Prepare worker arguments
    tasks = [(tree_file, reference_tips) for tree_file in replicate_trees]
    
    # Process replicates in parallel
    logger.info(f'Processing replicates with {cpus:,} CPUs.')
    final_split_counts = Counter()
    
    with Pool(processes=cpus) as pool:
        for split_counter in tqdm(
                pool.imap_unordered(process_tree, tasks),
                total=len(tasks),
                desc="Processing replicates",
                unit=" trees",
        ):
            final_split_counts.update(split_counter)
    
    # Apply support values to original tree
    logger.info('Applying bootstrap support values to original tree.')
    apply_support(orig_tree, final_split_counts, len(replicate_trees))

    # get bootstrap support for original rooting so we can use
    # this as the support for the unrooted trifurcation
    root_support = None
    if left_root_tips:
        left_root_tips = frozenset(left_root_tips)
        if left_root_tips in final_split_counts:
            root_support = final_split_counts[left_root_tips]
    
    if root_support is None and right_root_tips:
        right_root_tips = frozenset(right_root_tips)
        root_support = final_split_counts[right_root_tips]

    root_support_perc = None
    if root_support:
        root_support_perc = (root_support / len(replicate_trees)) * 100
    
    # Write output tree, forcing support value for the trifurcation
    # indicating an unrooted tree
    logger.info('Writing out tree with bootstrap support values.')
    newick_string = orig_tree.write(tree_format=0)

    if root_support_perc:
        newick_string = newick_string.replace(';', f'{root_support_perc};')
    
    with open(output_tree, 'w') as f:
        f.write(f'{newick_string}\n')


def bootstrap_alignment(msa, output_file, frac=1.0):
    """Bootstrap multiple sequence alignment.

    True bootstrapping requires subsampling an alignment,
    with replacement, to construct new alignments
    with the same length as the input alignment. The
    'frac' parameter allows shorter bootstrap alignments
    to be generated in order to reduce computational
    demands.

    Parameters
    ----------
    msa : d[seq_id] -> seq
      Full multiple sequence alignment.
    output_file : str
      File to write bootstrapped alignment.
    frac : float
      Fraction of alignment to subsample.
    """
    
    alignment_len = len(msa[list(msa.keys())[0]])
    sample_len = int(alignment_len * frac)
    cols = [random.randint(0, alignment_len - 1) for _ in range(sample_len)]

    with open(output_file, 'w') as fout:
        for seq_id, seq in msa.items():
            fout.write('>' + seq_id + '\n')
            for col in cols:
                fout.write(seq[col])
            fout.write('\n')

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

import sys
import logging
from collections import defaultdict

from biolib.common import is_float

import dendropy

"""
To do:
 1. There is a serious hack in taxonomic_consistency which should be resolved, but
     requires the viral and plasmid phylogenies to be taxonomically consistent.
"""


class Taxonomy(object):
    """Manipulation of Greengenes-style taxonomy files and strings.

    This class currently assumes a Greengenes-style taxonomy
    string with the following 7 taxonomic ranks:
      d__; c__; o__; f__; g__; s__

    Spaces after the semi-colons are optional.
    """

    rank_prefixes = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')
    rank_labels = ('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
    rank_index = {'d__': 0, 'p__': 1, 'c__': 2, 'o__': 3, 'f__': 4, 'g__': 5, 's__': 6}

    unclassified_rank = 'unclassified'

    unclassified_taxon = []
    for p in rank_prefixes:
        unclassified_taxon.append(p + unclassified_rank)
    unclassified_taxon = ';'.join(unclassified_taxon)

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()

    def taxa(self, tax_str):
        """Taxa specified by taxonomy string.

        Parameters
        ----------
        tax_str : str
            Greengenes-style taxonomy string.

        Returns
        -------
        list : [<domain>, <phylum>, ..., <species>]
            Rank order list of taxa.
        """

        taxa = [x.strip() for x in tax_str.split(';')]

        return taxa

    def taxa_at_ranks(self, tax_str):
        """Taxon at each taxonomic rank.

        Parameters
        ----------
        tax_str : str
            Greengenes-style taxonomy string.

        Returns
        -------
        dict : d[rank_label] -> taxon
            Taxon at each taxonomic rank.
        """

        taxa = self.taxa(tax_str)

        d = {}
        for rank, taxon in enumerate(taxa):
            d[Taxonomy.rank_labels[rank]] = taxon

    def check_full(self, tax_str):
        """Check if taxonomy string specifies all expected ranks.

        Parameters
        ----------
        tax_str : str
            Greengenes-style taxonomy string.

        Returns
        -------
        boolean
            True if string contains all expected ranks, else False.
        """

        taxa = [x.strip() for x in tax_str.split(';')]
        if len(taxa) != len(Taxonomy.rank_prefixes):
            self.logger.error('[Error] Taxonomy string contains too few ranks:')
            self.logger.error('[Error] %s' % str(taxa))
            return False

        for r, taxon in enumerate(taxa):
            if taxon[0:3] != Taxonomy.rank_prefixes[r]:
                self.logger.error('[Error] Taxon is not prefixed with the expected rank, %s.:' % Taxonomy.rank_prefixes[r])
                self.logger.error('[Error] %s' % str(taxa))
                return False

        return True

    def fill_missing_ranks(self, taxa):
        """Fill in any missing ranks in a taxonomy string.

        This function assumes the taxonomic ranks are
        in the proper rank order, but that some ranks may
        be missing.

        Parameters
        ----------
        list : [d__<taxon>, ..., s__<taxon>]
            List of taxa.

        Returns
        -------
        list
            List of taxa with all ranks.
        """

        new_taxa = []
        cur_taxa_index = 0
        for rank_prefix in Taxonomy.rank_prefixes:
            if len(taxa) > cur_taxa_index and taxa[cur_taxa_index][0:3] == rank_prefix:
                new_taxa.append(taxa[cur_taxa_index])
                cur_taxa_index += 1
            else:
                new_taxa.append(rank_prefix)

        return new_taxa

    def taxonomic_consistency(self, taxonomy, report_errors=True):
        """Determine taxonomically consistent classification for taxa at each rank.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        report_errors : boolean
            Flag indicating if errors should be written to screen.

        Returns
        -------
        dict : d[taxa] -> expected parent
            Expected parent taxon for taxa at all taxonomic ranks, or
            None if the taxonomy is inconsistent.
        """

        expected_parent = {}
        for genome_id, taxa in taxonomy.iteritems():
            if taxa[0] == 'd__Viruses' or '[P]' in taxa[0]:
                # *** This is a HACK. It would be far better to enforce
                # a taxonomically consistent taxonomy, but
                # the viral taxonomy at IMG is currently not consistent
                continue

            for r in xrange(1, len(taxa)):
                if taxa[r] == Taxonomy.rank_prefixes[r]:
                    break

                if taxa[r] in expected_parent:
                    if report_errors:
                        if taxa[r - 1] != expected_parent[taxa[r]]:
                            self.logger.error('[Error] Provided taxonomy is not taxonomically consistent.')
                            self.logger.error('[Error] Genome %s indicates the parent of %s is %s.' % (genome_id, taxa[r], taxa[r - 1]))
                            self.logger.error('[Error] The parent of this taxa was previously indicated as %s.' % (expected_parent[taxa[r]]))
                            return None

                expected_parent[taxa[r]] = taxa[r - 1]

        return expected_parent

    def extract_valid_species_name(self, taxon):
        """Try to extract a valid species name from a taxonomic label.

        A full species name should be  binomial and include a 'generic name' (genus) and
        a 'specific epithet' (species), i.e. Escherichia coli. This method
        assumes the two names should be separated by a space.

        Parameters
        ----------
        taxon : str
            Taxon label to process.

        Returns
        -------
        str
            Valid species name, or None.
        """

        if ' bacterium' in taxon.lower() or 'sp.' in taxon.lower():
            return None

        taxon = taxon.replace('s__', '')
        taxon = taxon.replace('Candidatus', '')
        taxon = taxon.replace('candidatus', '')

        if not taxon or taxon[0].islower():
            return None

        taxon_split = taxon.split(' ')
        if len(taxon_split) < 2:
            return None

        # sanity check
        taxon = 's__' + ' '.join(taxon_split[0:2])
        self.validate_species_name(taxon)

        return taxon

    def validate_species_name(self, species_name, require_full=True, require_prefix=True):
        """Validate species name.

        A full species name should be  binomial and include a 'generic name' (genus) and
        a 'specific epithet' (species), i.e. Escherichia coli. This method
        assumes the two names should be separated by a space.

        Parameters
        ----------
        species_name : str
            Species name to validate
        require_full : boolean
            Flag indicating if species name must include 'generic name and 'specific epithet'.
        require_prefix : boolean
            Flag indicating if name must start with the species prefix ('s__').

        Returns
        -------
        boolean
            True if species name is valid, otherwise False.
        str
            Reason for failing validation, otherwise None.
        """
        
        if species_name == 's__':
            return True, None
        
        # test for prefix
        if require_prefix:
            if not species_name.startswith('s__'):
                return False, 'name is missing the species prefix'
            
        # remove prefix before testing other properties
        test_name = species_name
        if test_name.startswith('s__'):
            test_name = test_name[3:]

        # test for full name
        if require_full:
            if 'candidatus' in test_name.lower():
                if len(test_name.split(' ')) <= 2:
                    return False, 'name appears to be missing the generic name'
            else:
                if len(test_name.split(' ')) <= 1:
                    return False, 'name appears to be missing the generic name'
	
	# check for tell-tale signs on invalid species names
	if " bacterium" in test_name.lower():
            return False, "name contains the word 'bacterium'"
	if " archaeon" in test_name.lower():
	    return False, "name contains the word 'archaeon'"
	if " archeaon" in test_name.lower():
	    return False, "name contains the word 'archeaon'"
	if "-like" in test_name.lower():
	    return False, "name contains '-like'"
        if " group " in test_name.lower():
	    return False, "name contains 'group'"
	if " symbiont" in test_name.lower():
	    return False, "name contains 'symbiont'"
	if " endosymbiont" in test_name.lower():
	    return False, "name contains 'endosymbiont'"
	if " taxon" in test_name.lower():
	    return False, "name contains 'taxon'"
	if " cluster" in test_name.lower():
	    return False, "name contains 'cluster'"
	if " of " in test_name.lower():
	    return False, "name contains 'of'"
	if test_name[0].islower():
	    return False, 'first letter of name is lowercase'
	if 'sp.' in test_name.lower():
	    return False, "name contains 'sp.'"

        return True, None

    def validate(self, taxonomy, check_prefixes, check_ranks, check_hierarchy, check_species, report_errors=True):
        """Check if taxonomy forms a strict hierarchy with all expected ranks.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        check_prefixes : boolean
            Flag indicating if prefix of taxon should be validated.
        check_ranks : boolean
            Flag indicating if the presence of all ranks should be validated.
        check_hierarchy : boolean
            Flag indicating if the taxonomic hierarchy should be validated.
        check_species : boolean
            Flag indicating if the taxonomic consistency of named species should be validated.
        report_errors : boolean
            Flag indicating if errors should be written to screen.

        Returns
        -------
        dict : d[taxon_id] -> taxonomy
            Taxa with invalid number of ranks.
        dict : d[taxon_id] -> [taxon, taxonomy]
            Taxa with invalid rank prefixes.
        dict: d[taxon_id] -> [species name, error message]
            Taxa with invalid species names.
        dict: d[child_taxon_id] -> two or more parent taxon ids
            Taxa with invalid hierarchies.
        """

        # check for incomplete taxonomy strings or unexpected rank prefixes
        invalid_ranks = {}
        invalid_prefixes = {}
        invalid_species_name = {}
        for taxon_id, taxa in taxonomy.iteritems():
            if check_ranks:
                if len(taxa) != len(Taxonomy.rank_prefixes):
                    invalid_ranks[taxon_id] = ';'.join(taxa)
                    continue

            if check_prefixes:
                for r, taxon in enumerate(taxa):
                    if taxon[0:3] != Taxonomy.rank_prefixes[r]:
                        invalid_prefixes[taxon_id] = [taxon, ';'.join(taxa)]
                        break

            if check_species:
                species_index = Taxonomy.rank_index['s__']
                if len(taxa) > species_index:
                    species_name = taxa[species_index]
                    valid, error_msg = self.validate_species_name(species_name, require_full=True, require_prefix=True)
                    if not valid:
                        invalid_species_name[taxon_id] = [species_name, error_msg]

        # check for inconsistencies in the taxonomic hierarchy
        invalid_hierarchies = defaultdict(set)
        if check_hierarchy:
            expected_parent = self.taxonomic_consistency(taxonomy, False)
            for taxon_id, taxa in taxonomy.iteritems():
                for r in xrange(1, len(taxa)):
                    if taxa[r] == Taxonomy.rank_prefixes[r]:
                        continue

                    if r == self.rank_index['s__'] and not check_species:
                        continue

                    if taxa[r - 1] != expected_parent[taxa[r]]:
                        invalid_hierarchies[taxa[r]].add(taxa[r - 1])
                        invalid_hierarchies[taxa[r]].add(expected_parent[taxa[r]])

        if report_errors:
            if len(invalid_ranks):
                print ''
                print 'Taxonomy contains too few ranks:'
                for taxon_id, taxa_str in invalid_ranks.iteritems():
                    print '%s\t%s' % (taxon_id, taxa_str)

            if len(invalid_prefixes):
                print ''
                print 'Taxonomy contains an invalid rank prefix:'
                for taxon_id, info in invalid_prefixes.iteritems():
                    print '%s\t%s\t%s' % (taxon_id, info[0], info[1])

            if len(invalid_species_name):
                print ''
                print 'Taxonomy contain invalid species names:'
                for taxon_id, info in invalid_species_name.iteritems():
                    print '%s\t%s\t%s' % (taxon_id, info[0], info[1])

            if len(invalid_hierarchies):
                print ''
                print 'Taxonomy contains taxa with multiple parents:'
                for child_taxon, parent_taxa in invalid_hierarchies.iteritems():
                    print '%s\t%s' % (child_taxon, ', '.join(parent_taxa))

        return invalid_ranks, invalid_prefixes, invalid_species_name, invalid_hierarchies

    def taxon_children(self, taxonomy):
        """Get children taxa for each taxonomic group.

        For species, this is a list of extant taxa.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[taxon] -> list of children taxa
            All children taxa for of each named taxonomic group.
        """

        taxon_children = defaultdict(set)
        for taxon_id, taxa in taxonomy.iteritems():
            for i, taxon in enumerate(taxa):
                if len(taxon) == 3:
                    continue  # just rank prefix

                if len(taxa) > i + 1 and len(taxa[i + 1]) != 3:
                    taxon_children[taxon].add(taxa[i + 1])

            if len(taxa) > self.rank_index['s__']:
                taxon = taxa[self.rank_index['s__']]
                if taxon != 's__':
                    taxon_children[taxon].add(taxon_id)

        return taxon_children

    def children(self, taxon, taxonomy):
        """Get children of taxon.

        For species, this is a list of extant taxa. For higher
        ranks, this is named groups and does not include the
        extant taxa.

        Parameters
        ----------
        taxon : str
            Named taxonomic group of interest.
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        set : {child1, child2, ..., childN}
            All children taxa for the named taxonomic group.
        """

        c = set()
        for taxon_id, taxa in taxonomy.iteritems():
            if taxon in taxa:

                if taxon.startswith('s__'):
                    c.add(taxon_id)
                else:
                    taxon_index = taxa.index(taxon)
                    for child in taxa[taxon_index + 1:]:
                        if len(child) > 3:  # not just an empty prefix
                            c.add(child)

        return c

    def parents(self, taxonomy):
        """Get parents for all taxa.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        d[taxon] -> list of parent taxa in rank order
            Parent taxa for each taxon.
        """

        p = defaultdict(list)
        for taxon_id, taxa in taxonomy.iteritems():
            for i, taxon in enumerate(taxa):
                p[taxon_id] = taxa
                if i != 0:
                    p[taxon] = taxa[0:i]

        return p

    def extant_taxa_for_rank(self, rank, taxonomy):
        """Get extant taxa for all named groups at the specified rank.

        Parameters
        ----------
        rank : str (e.g., class or order)
            Rank of interest
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[taxa] -> set of extant taxa
            Extant taxa for named groups at the specified rank.
        """

        assert(rank in Taxonomy.rank_labels)

        d = defaultdict(set)
        rank_index = Taxonomy.rank_labels.index(rank)
        for taxon_id, taxa in taxonomy.iteritems():
            if taxa[rank_index] != Taxonomy.rank_prefixes:
                d[taxa[rank_index]].add(taxon_id)

        return d

    def read_from_tree(self, tree):
        """Obtain the taxonomy for each extant taxa as specified by internal tree labels.

        Parameters
        ----------
        tree : str or dendropy.Tree
            Filename of newick tree or dendropy tree object.

        Returns
        -------
        dict : d[unique_id] -> [d__<taxon>, ..., s__<taxon>]
            Taxa indexed by unique ids.
        """

        if isinstance(tree, basestring):
            tree = dendropy.Tree.get_from_path(tree, schema='newick', rooting="force-rooted", preserve_underscores=True)

        taxonomy = {}
        for leaf in tree.leaf_node_iter():
            taxa = []

            node = leaf.parent_node
            while node:
                if node.label:
                    taxa_str = node.label
                    if ':' in taxa_str:
                        taxa_str = taxa_str.split(':')[1]

                    if not is_float(taxa_str):
                        # appears to be an internal label and not simply a support value
                        taxa = [x.strip() for x in taxa_str.split(';')] + taxa
                node = node.parent_node

            if len(taxa) > 7:
                self.logger.error('Invalid taxonomy string read from tree for taxon %s: %s' % (leaf.taxon.label, taxa))
                sys.exit(-1)

            # check if genus name should be appended to species label
            if len(taxa) == 7:
                genus = taxa[5][3:]
                species = taxa[6][3:]
                if genus not in species:
                    taxa[6] = 's__' + genus + ' ' + species

            taxa = self.fill_missing_ranks(taxa)
            taxonomy[leaf.taxon.label] = taxa

        return taxonomy

    def read(self, taxonomy_file):
        """Read Greengenes-style taxonomy file.

        Expected format is:
            <id>\t<taxonomy string>

        where the taxonomy string has the formats:
            d__; c__; o__; f__; g__; s__

        Parameters
        ----------
        taxonomy_file : str
            Greengenes-style taxonomy file.

        Returns
        -------
        dict : d[unique_id] -> [d__<taxon>, ..., s__<taxon>]
            Taxa indexed by unique ids.
        """

        d = {}
        for line in open(taxonomy_file):
            line_split = line.split('\t')
            unique_id = line_split[0]

            tax_str = line_split[1].rstrip()
            if tax_str[-1] == ';':
                # remove trailing semicolons which sometimes
                # appear in Greengenes-style taxonomy files
                tax_str = tax_str[0:-1]

            d[unique_id] = [x.strip() for x in tax_str.split(';')]

        return d

    def write(self, taxonomy, output_file):
        """Write Greengenes-style taxonomy file.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        output_file : str
            Name of output file.
        """

        fout = open(output_file, 'w')
        for genome_id, taxa in taxonomy.iteritems():
            fout.write(genome_id + '\t' + ';'.join(taxa) + '\n')
        fout.close()

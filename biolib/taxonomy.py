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


class InvalidTaxonomy(Exception):
    """Exception thrown for invalid taxonomy string."""
    pass
    

class Taxonomy(object):
    """Manipulation of Greengenes-style taxonomy files and strings.

    This class currently assumes a Greengenes-style taxonomy
    string with the following 7 taxonomic ranks:
      d__; c__; o__; f__; g__; s__
      
    Spaces after the semi-colons are optional.
    """

    def __init__(self):
        """Initialization."""
        self.rank_prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        self.rank_labels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        self.rank_index = {'d__': 0, 'p__': 1, 'c__': 2, 'o__': 3, 'f__': 4, 'g__': 5, 's__': 6]}
        
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

        for taxa
        
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
            d[self.rank_labels[i]] = taxon
            
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
        if len(taxa) != len(self.rank_prefixes):
            return False
            
        for rank, taxon in taxa:
            if taxon[0:3] != self.rank_prefixes[rank]:
                return False
                
        return True
        
    def fill_missing_ranks(self, tax_str):
        """Fill in any missing ranks in a taxonomy string.
        
        This function assumes the taxonomy string lists
        taxa in proper rank order, but that some ranks
        may be missing.
        
        Parameters
        ----------
        tax_str : str
            Greengenes-style taxonomy string.
        
        Returns
        -------
        str
            Taxonomy string with prefixes for all ranks.
        """
        
        taxa = [x.strip() for x in tax_str.split(';')]
        
        new_tax = []
        cur_taxa_index = 0
        for rank_prefix in self.rank_prefixes:
            if taxa[cur_taxa_index][0:3] == rank_prefix:
                cur_taxa_index.append(taxa[cur_taxa_index])
                cur_taxa_index += 1
            else:
                new_tax.append(rank_prefix)
                
        return ';'.join(new_tax)
        
    def read(self, taxonomy_file, validate=False):
        """Read Greengenes-style taxonomy file.
        
        Expected format is:
            <id>\t<taxonomy string>
        
        where the taxonomy string has the formats:
            d__; c__; o__; f__; g__; s__
        
        Parameters
        ----------
        taxonomy_file : str
            Greengenes-style taxonomy file.
        validate : boolean
            Check if all taxonomy strings are valid.
        
        Returns
        -------
        dict : d[unique_id] -> taxonomy string
            Taxonomy strings indexed by unique ids.
            
        Exceptions
        ----------
        InvalidTaxonomy 
            Thrown if validate is True and an invalid 
            taxonomy string is encountered.
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
            
            if validate and not self.check_full(tax_str):
                raise InvalidTaxonomy("Invalid taxonomy string: %s" % tax_str)
            
            d[unique_id] = tax_str
            
        return d
        
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

class Taxonomy(object):
    """Manipulation of Greengenes-style taxonomy strings.

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
        
    def split(self, taxa_str):
        """Split taxonomy string into ranks.
        
        Parameters
        ----------
        taxa_str : str
            Greengenes-style taxonomy string.
        
        Returns
        -------
        dict
            Taxon at each taxonomic rank.
        """
        
        taxa = [x.strip() for x in taxa_str.split(';')]
        
        d = {}
        for rank, taxon in enumerate(taxa):
            d[self.rank_labels[i]] = taxon
            
    def check_full(self, taxa_str):
        """Check if taxonomy string specifies all expected ranks.
        
        Parameters
        ----------
        taxa_str : str
            Greengenes-style taxonomy string.
        
        Returns
        -------
        boolean
            True if string contains all expected ranks, else False.
        """
        
        taxa = [x.strip() for x in taxa_str.split(';')]
        if len(taxa) != len(self.rank_prefixes):
            return False
            
        for rank, taxon in taxa:
            if taxon[0:3] != self.rank_prefixes[rank]:
                return False
                
        return True
        
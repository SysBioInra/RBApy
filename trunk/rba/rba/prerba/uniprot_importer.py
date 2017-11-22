"""Module defining UniprotImporter class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import urllib
try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen


class UniprotImporter(object):
    """
    Class retrieving uniprot data for specified organism.

    Attributes
    ----------
    data : ?
        Data downloaded from uniprot.

    """

    def __init__(self, organism_id):
        """
        Build object from uniprot organism identifier.

        Parameters
        ----------
        organism_id : str
            Information used to retrieve organism. It can be
            a uniprot identifier, a species name, etc.

        """
        # code adapted from
        # http://www.uniprot.org/help/programmatic_access
        url = 'http://www.uniprot.org/uniprot/'
        params = {
            'format': 'tab',
            'query': 'organism:' + organism_id,
            'columns': url_columns()
        }
        url_data = urllib.urlencode(params)
        request = Request(url, url_data)
        contact = "anne.goelzer@inra.fr"
        request.add_header('User-Agent', 'Python %s' % contact)
        response = urlopen(request)
        self.data = response.read()


def url_columns():
    """Build the url part that specifies which columns we want."""
    # Information needed to build url was retrieved on
    # http://www.uniprot.org/help/uniprotkb_column_names

    ######################
    # Names and taxonomy #
    ######################
    # **Entry**, **Entry name**,
    # **Gene names**, **Protein names**, **Organism**, **Organism ID**
    cols = ['id', 'entry name', 'genes', 'protein names', 'organism',
            'organism-id']

    #############
    # Sequences #
    #############
    # **Length**, **Mass**, **Sequence**
    cols += ['length', 'mass', 'sequence']

    ############
    # Function #
    ############
    # **EC number**, **Catalytic activity**, **Cofactor**,
    # **Enzyme regulation**, **Function [CC]**, **Pathway**
    # **Temperature dependence**, **pH dependence**
    # **Metal binding**, **Nucleotide binding**
    comments = ['catalytic activity', 'cofactor', 'enzyme regulation',
                'function', 'pathway', 'temperature dependence',
                'ph dependence']
    features = ['metal binding', 'np binding']
    cols += (['ec'] + reformat('comment', comments)
             + reformat('feature', features))

    #################
    # Miscellaneous #
    #################
    # **Features**, **Caution**, **Keywords**
    cols += ['features'] + reformat('comment', ['caution'] + ['keywords'])

    ###############
    # Interaction #
    ###############
    # **Subunit structure [CC]**
    cols += reformat('comment', ['subunit'])

    ##############
    # Expression #
    ##############
    # **Tissue specificity**
    cols += reformat('comment', ['tissue specificity'])

    ########################
    # Subcellular location #
    ########################
    # **Subcellular location [CC]**
    cols += reformat('comment', ['subcellular location'])

    return ','.join(cols)


def reformat(keyword, fields):
    """
    Reformat field name to url format using specific keyword.

    Example:
        reformat('comment', ['a','b']) returns
        ['comment(A)', 'comment(B)']
    """
    return ['{}({})'.format(keyword, f.upper()) for f in fields]

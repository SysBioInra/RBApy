"""Module defining UniprotImporter class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import os.path
import sys
import requests
#try:
#    from urllib.request import urlopen
#    from urllib.parse import urlencode
#except ImportError:
#    from urllib2 import urlopen
#    from urllib import urlencode


def create_uniprot_if_absent(input_file, organism_id):
    if not os.path.isfile(input_file):
        print('Could not find UniProt file. Downloading most recent'
              ' version ...', end='')
        sys.stdout.flush()
        raw_data = UniprotImporter(organism_id).data
        if len(raw_data) == 0:
            raise UserWarning('Invalid organism, could not retrieve '
                              'UniProt data.')
        with open(input_file, 'wb') as f:
            f.write(raw_data)
        print(' done')


class UniprotImporter(object):
    """
    Class retrieving UniProt data for specified organism.

    Attributes
    ----------
    data : ?
        Data downloaded from UniProt.

    """

    def __init__(self, organism_id):
        """
        Build object from UniProt organism identifier.

        Parameters
        ----------
        organism_id : str
            Information used to retrieve organism. It can be
            a UniProt identifier, a species name, etc.

        """
        # code adapted from
        # http://www.UniProt.org/help/programmatic_access
        #url = 'http://www.uniprot.org/uniprot/'
        #params = {
        #    'format': 'tab',
        #    'query': 'organism:' + organism_id,
        #    'columns': url_columns()
        #}
        #url_data = urlencode(params)
        #response = urlopen(url + '?' + url_data)
        #self.data = response.read()

        fields_to_look_up=['accession','id','gene_names','annotation_score','protein_name','organism_name','organism_id','reviewed','length','mass','sequence','ec',
        'cc_catalytic_activity','cc_cofactor','cc_subunit','cc_subcellular_location','cc_tissue_specificity','cc_pathway','cc_function',
        'temp_dependence','ph_dependence','cc_caution','cc_interaction','cc_ptm','ft_signal','ft_transit','ft_propep','ft_lipid','ft_carbohyd','ft_disulfid','ft_intramem','ft_transmem','redox_potential']
        #testwise removed:
        # fields unable to be retreived: 'ft_metal' and 'ft_np_bind'

        payload = {'query':'organism_id:{}'.format(organism_id),
                   'format':'tsv',
                   'fields':fields_to_look_up}
        url = 'https://rest.uniprot.org/uniprotkb/stream'
        response = requests.get(url, params=payload)
        self.data = response.content

# obsolete?
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

    ######################
    # Annotation quality #
    ######################
    # **Reviewed**, **Annotation score**
    cols += ['reviewed', 'annotation score']


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

# obsolete?
def reformat(keyword, fields):
    """
    Reformat field name to url format using specific keyword.

    Example:
        reformat('comment', ['a','b']) returns
        ['comment(A)', 'comment(B)']
    """
    return ['{}({})'.format(keyword, f.upper()) for f in fields]

"""Module defining UniprotImporter class."""

# python 2/3 compatibility
from __future__ import division, print_function, absolute_import

# global imports
import os.path
import sys
import requests


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
        #fields_to_look_up=['accession','id','gene_names','annotation_score','protein_name','organism_name','organism_id','reviewed','length','mass','sequence','ec',
        #'cc_catalytic_activity','cc_cofactor','cc_subunit','cc_subcellular_location','cc_tissue_specificity','cc_pathway','cc_function',
        #'temp_dependence','feature_count','ph_dependence','cc_caution','cc_interaction','cc_ptm','ft_signal','ft_transit','ft_propep','ft_lipid','ft_carbohyd','ft_disulfid','ft_intramem','ft_transmem','redox_potential']
        #fields unable to be retreived: 'ft_metal' and 'ft_np_bind'

        payload = {'query':'organism_id:{}'.format(organism_id),
                   'format':'tsv',
                   'fields':url_columns()}
        url = 'https://rest.uniprot.org/uniprotkb/stream'
        response = requests.get(url, params=payload)
        self.data = response.content

def url_columns():
    """Build the url part that specifies which columns we want."""
    # Information needed to build url was retrieved on
    # http://www.uniprot.org/help/uniprotkb_column_names

    ######################
    # Names and taxonomy #
    ######################
    # **Entry**, **Entry Name**,
    # **Gene Names**, **Protein names**, **Organism**, **Organism (ID)**
    cols= ['accession','id','gene_names','protein_name','organism_name','organism_id']

    ######################
    # Annotation quality #
    ######################
    # **Reviewed**, **Annotation**
    cols += ['reviewed', 'annotation_score']

    #############
    # Sequences #
    #############
    # **Length**, **Mass**, **Sequence**
    cols += ['length', 'mass', 'sequence']

    ############
    # Function #
    ############
    # **Catalytic activity** , **Pathway** ,
    # **Function [CC]** , **EC number** ,
    # **Cofactor** , **Temperature dependence** ,
    # **pH dependence** , **Redox potential**
    cols += ['cc_catalytic_activity','cc_pathway','cc_function','ec','cc_cofactor','temp_dependence','ph_dependence','redox_potential']
    # **Metal binding**, **Nucleotide binding**
    # fields unable to be retreived: 'ft_metal' and 'ft_np_bind'

    ###############
    # Interaction #
    ###############
    # **Subunit structure** , **Interacts with**
    cols += ['cc_subunit','cc_interaction']

    ########################
    # Subcellular location #
    ########################
    # **Subcellular location [CC]** , **Intramembrane** , **Transmembrane**
    cols += ['cc_subcellular_location','ft_intramem','ft_transmem']

    ##############
    # Expression #
    ##############
    # **Tissue specificity**
    cols += ['cc_tissue_specificity']

    ##############
    # Protein processing #
    ##############
    # **Post-translational modification** , **Lipidation**
    # **Glycosylation** , **Disulfide bond**
    cols += ['cc_ptm','ft_lipid','ft_carbohyd','ft_disulfid']

    ##############
    # Specific peptides #
    ##############
    # **Signal peptide** , **Transit peptide** , **Propeptide**
    cols += ['ft_signal','ft_transit','ft_propep']

    #################
    # Miscellaneous #
    #################
    #  **Caution**,**Features**, **Keywords**
    cols += ['cc_caution','feature_count','keyword']

    return cols

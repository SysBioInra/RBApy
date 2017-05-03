def postprocess_sbml(sbml_data):
    """
    Postprocess protein names in SBML data.

    :param sbml_data: Filtered SBML data.
    :type sbml_data: SBMLFilter instance.
    """
    # remove fake enzyme components
    # diffusion
    sbml_data.remove_protein_matching(r'^d\w+$')
    # spontaneous
    sbml_data.remove_protein_matching(r'^s\w+$')
    # RNAs ??? <-- check this
    sbml_data.remove_protein_matching(r'^RS0\w+$')
    # replace unknown proteins with average proteins
    # unknown excretion proteins
    sbml_data.replace_protein_matching(r'^e\w+$', 'AverageEnzymaticProtein')
    # unknown proteins ???
    sbml_data.replace_protein_matching(r'NoAssignment$', 'AverageEnzymaticProtein')
    # no_assignment_reaction = my_sbml.reactions_associated_to_protein('NoAssignment')

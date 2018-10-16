"""

Functions for performing compatibility checks and conversions for MD programs.

"""

def check_term_compatibility(dl, program_term_data):
    """
    Checks that term parameters in datalayer are compatible with the output program

    Parameters
    -------------
    dl : eex datalayer
        The EEX datalayer object containing information for the system being translated

    program_metadata: dict
        Term metadata for program being converted to


    Return
    ---------------
    dl : eex datalayer
        EEX datalayer which is compatible with program metadata
    """



    # Check that keywords are compatible
    for order, term_data in program_term_data.items():

        stored_terms = dl.get_terms(order)
        stored_keywords = stored_terms.keys()
        compatible_keywords = term_data.keys()

        #print(compatible_keywords, stored_keywords)




    # Perform conversion if necessary
    # Return modified dl

    return dl

def check_atom_metadata():
    """
    Check that all necessary atom metadata is present.
    :return:
    """


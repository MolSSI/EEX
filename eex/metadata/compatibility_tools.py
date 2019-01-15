"""

Functions for performing compatibility checks and conversions for MD programs.

"""

import eex
import numpy as np

def check_term_compatibility(dl, program_term_data):
    """
    Checks that term parameters in datalayer are compatible with the output program. Performs conversion if necessary
    and possible.

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

    # Check that keywords are compatible. This loops through the term data of the program to convert to
    for order, term_data in program_term_data.items():

        # Get allowed keywords for conversion (this is what is in the metadata for a program)
        allowed_keywords = list(term_data.keys())
        allowed_canonical_keywords = []

        # Get canonical forms for what is "allowed"
        for allowed_keyword in allowed_keywords:
            # Need utility function for this instead of this mess. TODO
            term_md = eex.metadata.get_term_metadata(order, "forms", allowed_keyword)
            allowed_canonical_keywords.append(term_md["canonical_form"])

        # Get stored keywords and uids
        stored_terms = dl.list_term_parameters(order)
        uids = list(stored_terms.keys())

        stored_keywords = []
        [stored_keywords.append(stored_terms[uids[x]][0]) for x in uids]

        stored_canonical_keywords = []

        # Get canonical forms for what is stored.
        for stored_keyword in stored_keywords:
            term_md = eex.metadata.get_term_metadata(order, "forms", stored_keyword)
            stored_canonical_keywords.append(term_md["canonical_form"])

        # Get difference between what is stored in datalayer and what is allowed. Diff will give you what is in
        # the datalayer. "Canonical" forms and those that are stored are both compared here.
        diff_canonical = np.setdiff1d(stored_canonical_keywords, allowed_canonical_keywords)
        diff = np.setdiff1d(stored_keywords, allowed_keywords)

        # If diff is found between canonical, not compatible. Otherwise try for conversion.
        if diff_canonical.size > 0:
            raise TypeError("Functional forms, %s, found in datalayer are not compatible with required forms %s" % (
            diff_canonical, allowed_keywords))
        elif diff.size == 0:
            # Then no conversion is necessary
            pass
        else:
            # The conversion here is complicated/confusing and can likely be improved to be much clearer with some utility functions.

            for uid, term in list(stored_terms.items()):
                # Need a utility function for this. TODO
                term_parameters = eex.metadata.get_term_metadata(order, "forms", term[0])["parameters"]
                convert_parameters = {k: v for k,v in zip(term_parameters, term[1:])}

                # Call conversion function. TODO - this conversion may not work (allowed_keywords[0])
                print("Converting",  term[0], "to", allowed_keywords[0])
                new_form = eex.form_converters.convert_form(order, convert_parameters, term[0], allowed_keywords[0])

                # Remove old form

                # Get indices for terms to remove
                all_terms = dl.get_terms(order)
                remove_terms = all_terms[all_terms["term_index"] == uid]
                atom_columns = [x for x in remove_terms.columns if "atom" in x]

                dl.remove_terms(order, index=remove_terms.index, propagate=False)

                # Remove this term parameter from datalayer. TODO - problem here for multiple removals!!
                dl.remove_term_parameter(order, uid)

                # Add new form - Two things must be done here:
                # First, the converted functional form must be added using dl.add_term_parameter
                # Then, terms themselves must be re-added using dl.add term

                # Build dataframe
                # 1. Get atoms involed in term from remvoe terms - 'base_atoms'
                # 2. Build df with atoms and term_index

                base_atoms = remove_terms[atom_columns]
                new_form_keys = new_form.keys()

                # Loop through new terms to add to dl
                for val in range(len(list(new_form.values())[0])):
                    to_add = {}

                    for key in new_form_keys:
                        to_add[key] = new_form[key][val]

                    uid = dl.add_term_parameter(order, allowed_keywords[0], to_add)
                    base_atoms['term_index'] = uid

                    dl.add_terms(order, base_atoms)

    # Return modified dl

    return dl

def check_atom_metadata():
    """
    Check that all necessary atom metadata is present.
    :return:
    """


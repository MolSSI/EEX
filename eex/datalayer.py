"""
Contains the DataLayer class (name in progress) which takes and reads various peices of data
"""

import os
import pandas as pd
import numpy as np
import collections
import copy

import eex
from . import filelayer
from . import metadata
from . import units
from . import utility

APC_DICT = metadata.atom_property_to_column


class DataLayer(object):
    def __init__(self, name, store_location=None, save_data=False, backend="Memory"):
        """
        Initializes the DataLayer class

        Parameters
        ----------
        name : str
            The name of the energy expression stored
        store_location : {None, str}, optional
            The location to store the temporary data during the translation. Defaults to the current working directory.
        save_data : {False, True}, optional
            Decides whether to delete the store data upon destruction of the DataLayer object.
        backend : {"HDF5", "memory"}, optional
            Storage backend for the energy expression.
        """

        # Set the state
        self.name = name

        # Build the store
        self.store_location = store_location
        if self.store_location is None:
            self.store_location = os.getcwd()

        if backend.upper() == "HDF5":
            self.store = filelayer.HDFStore(self.name, self.store_location, save_data)
        elif backend.upper() == "MEMORY":
            self.store = filelayer.MemoryStore(self.name, self.store_location, save_data)
        else:
            raise KeyError("DataLayer: Backend of type '%s' not recognized." % backend)

        # Setup empty parameters dictionary
        self._functional_forms = {2: {}, 3: {}, 4: {}}
        self._terms = {2: {}, 3: {}, 4: {}}
        self._atom_metadata = {}
        self._atom_sets = set()

### Generic helper close/save/list/etc functions

    def call_by_string(self, *args, **kwargs):
        """
        Adds the ability to call DL function by their string name.
        """

        if args[0] == "NYI":
            return False

        try:
            function = getattr(self, args[0])
        except:
            raise KeyError("DataLayer:call_by_string: does not have method %s." % args[0])

        return function(*args[1:], **kwargs)

    def close(self):
        """
        Closes the DL object
        """
        self.store.close()

    def list_tables(self):
        """
        Lists tables loaded into the store.
        """
        return [x for x in self.store.list_tables() if not x.startswith("other_")]

    def list_other_tables(self):
        """
        Lists "other" tables loaded into the store.
        """
        return [x.replace("other_", "") for x in self.store.list_tables() if x.startswith("other_")]

### Atom functions

    def _check_atoms_dict(self, property_name):
        if property_name in list(self._atom_metadata):
            return False

        field_data = metadata.atom_metadata[property_name]
        self._atom_metadata[property_name] = {"uvals": {}, "inv_uvals": {}}

        return True

    def _find_unqiue_atom_values(self, df, property_name):
        """
        Hashes the input parameters to build in internal index of unique values.
        """

        field_data = metadata.atom_metadata[property_name]
        param_dict = self._atom_metadata[property_name]

        cols = field_data["required_columns"]

        if field_data["dtype"] == float:
            df = df[cols].round(field_data["tol"])

        ret_df = pd.DataFrame(index=df.index)
        ret_df[property_name] = 0

        # For each unique value
        for gb_idx, udf in df.groupby(cols):

            # Update dictionary if necessary
            if gb_idx not in list(param_dict["uvals"]):

                # Bidirectional dictionary
                new_key = utility.find_lowest_hole(list(param_dict["inv_uvals"]))
                param_dict["uvals"][gb_idx] = new_key
                param_dict["inv_uvals"][new_key] = gb_idx

            # Grab the unique and set
            uidx = param_dict["uvals"][gb_idx]
            ret_df.loc[udf.index, property_name] = uidx

        return ret_df

    def _build_atom_values(self, df, property_name):
        """
        Expands the unique parameters using the built in property_name dictionary.
        """
        field_data = metadata.atom_metadata[property_name]
        param_dict = self._atom_metadata[property_name]

        cols = field_data["required_columns"]

        ret_df = pd.DataFrame(index=df.index)
        ret_df[property_name] = 0.0

        for gb_idx, udf in df.groupby(cols):
            ret_df.loc[udf.index, property_name] = param_dict["inv_uvals"][gb_idx]

        return ret_df

    def _store_atom_table(self, table_name, df, property_name, by_value, utype):
        """
        Internal way to store atom tables
        """

        self._check_atoms_dict(property_name)
        field_data = metadata.atom_metadata[property_name]

        # Figure out unit scaling factors
        if by_value and (field_data["units"] is not None) and (utype is not None):
            scale_factor = units.conversion_factor(utype, field_data["utype"])
            df = df[field_data["required_columns"]] * scale_factor

        if by_value and not (metadata.atom_metadata[property_name]["unique"]):
            tmp_df = self._find_unqiue_atom_values(df, property_name)
        else:
            tmp_df = df[APC_DICT[property_name]]

        return self.store.add_table(table_name, tmp_df)

    def _get_atom_table(self, table_name, property_name, by_value, utype):

        tmp = self.store.read_table(table_name)
        if by_value and not (metadata.atom_metadata[property_name]["unique"]):
            tmp = self._build_atom_values(tmp, property_name)

        # Figure out unit scaling factors
        field_data = metadata.atom_metadata[property_name]
        if by_value and (field_data["units"] is not None) and (utype is not None):
            scale_factor = units.conversion_factor(field_data["utype"], utype)
            tmp[field_data["required_columns"]] *= scale_factor

        return tmp

    def add_atom_parameters(self, property_name, value, uid=None, utype=None):

        property_name = property_name.lower()
        self._check_atoms_dict(property_name)
        param_dict = self._atom_metadata[property_name]
        field_data = metadata.atom_metadata[property_name]

        if (utype is not None) and (field_data["units"] is not None):
            value = value * units.conversion_factor(field_data["utype"], utype)

        if field_data["dtype"] == float:
            value = round(value, field_data["tol"])

        # Check if we have this key
        found_key = None
        for k, v in param_dict["inv_uvals"].items():
            if v == value:
                found_key = k

        # Brand new uid, return next in sequence
        if uid is None:
            if found_key is not None:
                return found_key

            new_key = utility.find_lowest_hole(list(param_dict["inv_uvals"]))
            param_dict["uvals"][value] = new_key
            param_dict["inv_uvals"][new_key] = value
            return new_key

        # We have a uid
        else:
            # Fine if it matches internally, otherwise throw
            if (found_key is not None):
                if (found_key == uid):
                    return found_key
                else:
                    raise KeyError(
                        "DataLayer:add_atom_parameters: Tried to add value %s, but found in uid (%d) and current keys (%d)"
                        % (value, uid, found_key))

            param_dict["inv_uvals"][uid] = value
            param_dict["uvals"][value] = uid
            return uid

    def list_atom_properties(self):
        """
        Lists all the valid atom properties that have been added.

        """
        return list(self._atom_sets)

    def add_atoms(self, atom_df, property_name=None, by_value=False, utype=None):
        """
        Adds atom information to the DataLayer object.

        Parameters
        ----------
        atom_df : {DataFrame, list, tuple}
            The atom data to add to the object.
        property_name : {list, str}, optional
            The atom property that is added, only necessary if a list is passed in.
        by_value : bool
            If data is passed by_value the DL automatically hashes the parameters to unique components.
        utype : {dict, pint.Unit}
            The unit type of the atom_df.

        Returns
        -------
        return : bool
            If the add was successful or not.

        Example
        -------
        dl = DataLayer("test")

        # Add the XYZ information for five random atoms
        tmp_df = pd.DataFrame(np.random.rand(5, 3), columns=["X", "Y", "Z"])
        tmp_df["atom_index"] = np.arange(5)
        dl.add_atom(tmp_df)
        """

        # Our index name
        index = "atom_index"

        # Validate DataFrame
        if not isinstance(atom_df, pd.DataFrame):
            raise KeyError("DataLayer:add_atoms: Data type '%s' not understood." % type(atom_df))

        if index in atom_df.columns:
            atom_df = atom_df.set_index(index, drop=True)

        if atom_df.index.name != index:
            raise KeyError("DataLayer:add_atoms: DF index must be the `atom_index` not '%s'." % atom_df.index.name)

        if utype is None:
            utype = {}
        elif isinstance(utype, dict):
            utype = {k.lower(): v for k, v in utype.items()}
        else:
            raise TypeError("utype type not understood")

        # Try to add all possible properties
        set_cols = set(atom_df.columns)
        found_one = False
        for k, v in APC_DICT.items():
            # Check if v is in the set_cols (set logic)
            if set(v) <= set_cols:
                uval = None
                if k in utype:
                    uval = utype[k]
                self._store_atom_table(k, atom_df, k, by_value, uval)
                found_one = True
                # Update what we have
                self._atom_sets |= set([k])
        if not found_one:
            raise Exception("DataLayer:add_atom: No data was added as no key was matched from input columns:\n%s" %
                            (" " * 11 + str(atom_df.columns)))

        return True

    def get_atoms(self, properties, by_value=False, utype=None):
        """
        Obtains atom information to the DataLayer object.

        Parameters
        ----------
        properties : {list, str}
            The properties to obtain for the atom data.
        by_value : bool
            If true returns the property by value, otherwise returns by index.

        Returns
        -------
        return : pd.DataFrame
            Returns a DataFrame containing the atom property information
            If the add was successful or not.

        """

        valid_properties = list(metadata.atom_property_to_column)

        if properties is None:
            properties = self.list_atom_properties()

        if not isinstance(properties, (tuple, list)):
            properties = [properties]

        # Make sure they are lower case
        properties = [x.lower() for x in properties]

        if not set(properties) <= set(list(valid_properties)):
            invalid_props = set(properties) - set(list(valid_properties))
            raise KeyError("DataLayer:add_atoms: Property name(s) '%s' not recognized." % str(list(invalid_props)))

        if utype is None:
            utype = {}
        elif isinstance(utype, dict):
            utype = {k.lower(): v for k, v in utype.items()}
        else:
            raise TypeError("utype type not understood")

        df_data = []
        for prop in properties:
            uval = None
            if prop in utype:
                uval = utype[prop]
            tmp = self._get_atom_table(prop, prop, by_value, uval)
            df_data.append(tmp)

        return pd.concat(df_data, axis=1)

### Term functions

    def add_parameters(self, order, term_name, term_parameters, uid=None, utype=None):
        """
        Adds the parameters of a registered functional form to the Datalayer object

        Parameters
        ----------
        order : int
            The order of the functional form (2, 3, 4, ...)
        term_name : str
            The name of the functional form you are adding.
        term_parameters : {list, tuple, dict}
            The parameters to the functional form you are adding. If a list or tuple the order matches the order supplied
            in the functional form. Otherwise the dictionary matches functional form parameter names.
        uid : int, optional
            The uid to assign to this parameterized term.
        utype : list of Pint units, options
            Custom units for this particular addition, otherwise uses the default units in the registered functional form.

        Examples
        --------

        assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0])
        assert 0 == dl.add_parameters(2, "harmonic", [4.0, 5.0])
        assert 1 == dl.add_parameters(2, "harmonic", [4.0, 6.0])

        """

        user_order = order
        order = metadata.sanitize_term_order_name(order)

        # Make sure we know what this is
        try:
            term_md = metadata.get_term_metadata(order, "forms", term_name)
        except KeyError:
            raise KeyError("DataLayer:add_parameters: Did not understand term '%d, %s'." % (order, term_name))

        # Validate and converate data as needed
        params = metadata.validate_term_dict(term_name, term_md, term_parameters, utype=utype)

        # First we check if we already have it
        found_key = None
        for k, v in self._terms[order].items():
            if (v[0] == term_name) and np.allclose(v[1:], params):
                found_key = k
                break

        # Figure out what actually to do
        if uid is None:

            if found_key is not None:
                return found_key

            # We have a new parameter! Find the lowest number that we can add it at.
            if not len(self._terms[order]):
                new_key = 0
            else:
                new_key = utility.find_lowest_hole(self._terms[order])

            params.insert(0, term_name)
            self._terms[order][new_key] = params

            return new_key

        else:

            if not isinstance(uid, int):
                raise TypeError(
                    "DataLayer:register_functional_forms: uid keyword must be of type int, found type '%s'." %
                    type(uid))

            # If we exist this could get dangerous
            if uid in self._terms[order]:
                old_param = self._terms[order][uid]
                match = (old_param[0] == term_name) and np.allclose(old_param[1:], params)
                if not match:
                    raise KeyError(
                        "DataLayer:register_functional_forms: uid already exists, but does not much current parameters."
                    )
                else:
                    return uid

            else:
                params.insert(0, term_name)
                self._terms[order][uid] = params

                return uid

    def get_parameters(self, order, uid, utype=None):

        order = metadata.sanitize_term_order_name(order)

        if (uid not in self._terms[order]):
            raise KeyError("DataLayer:get_parameters: Did not find term '%d %d" % (order, uid))

        # Stored as [term_name, parameters...]
        data = self._terms[order][uid]

        term_md = metadata.get_term_metadata(order, "forms", data[0])

        # Zip up the parameters
        parameters = {k: v for k, v in zip(term_md["parameters"], data[1:])}

        # Were done
        if utype is None:
            return (data[0], parameters)

        # Need to convert
        if isinstance(utype, (list, tuple)):
            if len(utype) != len(term_md["parameters"]):
                raise KeyError("DataLayer:get_parameters: length of utype should match the length of parameters.")
            utype = {k: v for k, v in zip(term_md["parameters"], utype)}

        if not isinstance(utype, dict):
            raise TypeError("DataLayer:get_parameters: Input utype '%s' is not understood." % str(type(utype)))

        for key in term_md["parameters"]:
            parameters[key] *= units.conversion_factor(term_md["utype"][key], utype[key])

        return (data[0], parameters)

    def list_parameter_uids(self, order=None):

        # Return everything
        if order is None:
            ret = {}
            for k, v in self._terms.items():
                ret[k] = list(v)
            return ret

        # Just return a specific order
        order = metadata.sanitize_term_order_name(order)

        return list(self._terms[order])

    def add_terms(self, order, df):

        order = metadata.sanitize_term_order_name(order)
        if order not in list(self._functional_forms):
            raise KeyError("DataLayer:add_terms: Did not understand order key '%s'." % str(order))

        req_cols = metadata.get_term_metadata(order, "index_columns")

        not_found = set(req_cols) - set(df.columns)
        if not_found:
            raise KeyError("DataLayer:add_terms: Missing required columns '%s' for order %d" % (str(not_found), order))

        if "term_index" in df.columns:
            df = df[req_cols + ["term_index"]]
        else:
            raise Exception("NYI: Add terms by *not* term_index")
        self.store.add_table("term" + str(order), df)

    def get_terms(self, order):
        order = metadata.sanitize_term_order_name(order)
        if order not in list(self._functional_forms):
            raise KeyError("DataLayer:add_terms: Did not understand order key '%s'." % str(order))

        return self.store.read_table("term" + str(order))

    def add_bonds(self, bonds):
        """
        Adds bond using a index notation.

        Parameters
        ----------
        bonds : pd.DataFrame
            Adds a DataFrame containing the bond information by index
            Required columns: ["bond_index", "atom1_index", "atom2_index", "bond_type"]

        Returns
        -------
        return : bool
            Returns a boolean value if the operations was successful or not
        """

        self.add_terms("bonds", bonds)

        return True

    def get_bonds(self):

        return self.get_terms("bonds")

    def add_angles(self, angles):

        self.add_terms("angles", angles)

    def get_angles(self):

        return self.get_terms("angles")

    def add_dihedrals(self, dihedrals):

        self.add_terms("dihedrals", dihedrals)

    def get_dihedrals(self):

        return self.get_terms("dihedrals")

### Other quantities

    def add_other(self, key, df):
        """
        Adds arbitrary data to the DataLayer object. This data is effectively private and will not be used by any part
        of EEX.

        Parameters
        ----------
        key : str
            The key to store the data under.
        df : pd.DataFrame
            Adds a DataFrame containing data

        Returns
        -------
        return : bool
            Returns a boolean value if the operations was successful or not
        """

        key = "other_" + key
        self.store.add_table(key, df)

        return True

    def get_other(self, key):
        """
        Obtains other information from the DataLayer object.
        """

        if not isinstance(key, (tuple, list)):
            key = [key]

        tmp_data = []
        for k in key:
            k = "other_" + k
            tmp_data.append(self.store.read_table(k))

        return pd.concat(tmp_data, axis=1)

"""
Contains the DataLayer class (name in progress) which takes and reads various pieces of data
"""

import copy
import json
import os

import numpy as np
import pandas as pd

from . import energy_eval
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

        self.store = filelayer.build_store(backend, self.name, self.store_location, save_data)

        # Setup empty term holder
        self._terms = {order: {} for order in [2, 3, 4]}
        self._term_count = {order: {"total": 0} for order in [2, 3, 4]}

        # Setup atom holders
        self._atom_metadata = {}
        self._atom_counts = {}

        for k, v in metadata.atom_metadata.items():
            if not v["unique"]:
                self._atom_metadata[k] = {"uvals": {}, "inv_uvals": {}}
        self._atom_counts = {k: 0 for k in list(metadata.atom_metadata)}

        # Any remaining metadata
        self._box_size = {}

### Generic helper close/save/list/etc functions

    def call_by_string(self, *args, **kwargs):
        """
        Adds the ability to call DL function by their string name.
        """

        if args[0] == "NYI":
            return False

        try:
            func = getattr(self, args[0])
        except AttributeError:
            raise AttributeError("DataLayer:call_by_string: does not have method %s." % args[0])

        return func(*args[1:], **kwargs)

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

    def set_box_size(self, bsize, utype=None):
        """
        Sets the overall size of the box for the datalayer
        """

        cf = 1.0
        if utype is not None:
            internal_length = units.convert_contexts("[length]")
            cf = units.conversion_factor(utype, internal_length)

        for key in ["x", "y", "z"]:
            if key.lower() in bsize:
                tmp = bsize[key]
            elif key.upper() in bsize:
                tmp = bsize[key].lower()
            else:
                raise KeyError("Could not find key '%s'." % key)

            if len(tmp) != 2:
                raise IndexError("bsize['%s'] length does not equal 2" % key)

            self._box_size[key] = (bsize[key][0] * cf, bsize[key][1] * cf)

    def get_box_size(self, utype=None):
        """
        Gets the overall size of the box for the datalayer
        """
        ret = copy.deepcopy(self._box_size)

        if utype is not None:
            internal_length = units.convert_contexts("[length]")
            cf = units.conversion_factor(internal_length, utype)
            for k, v in ret.items():
                ret[k] = (v[0] * cf, v[1] * cf)

            return ret

        else:
            return ret

    def evaluate(self):
        """
        Evaluate the current state of the energy expression.
        """
        return energy_eval.evaluate_energy_expression(self)

### Atom functions

    def _check_atoms_dict(self, property_name):
        """
        Builds the correct data struct in the atom metadata
        """
        property_name = property_name.lower()
        if property_name not in metadata.atom_metadata:
            raise Exception("DataLayer: Atom property %s is not valid." % property_name)

        return property_name

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

        field_data = metadata.atom_metadata[property_name]

        # Figure out unit scaling factors
        if by_value and (field_data["units"] is not None) and (utype is not None):
            scale_factor = units.conversion_factor(utype, field_data["utype"])
            df = df[field_data["required_columns"]] * scale_factor

        # Make sure our ints are ints and not accidentally floats
        if field_data["dtype"] == int:
            df = df[field_data["required_columns"]].astype(int, copy=True)

        # Handle the unique or filter data
        if by_value and not (metadata.atom_metadata[property_name]["unique"]):
            tmp_df = self._find_unqiue_atom_values(df, property_name)
        else:
            tmp_df = df[APC_DICT[property_name]]

        self._atom_counts[property_name] += tmp_df.shape[0]

        return self.store.add_table(table_name, tmp_df)

    def _get_atom_table(self, table_name, property_name, by_value, utype):

        tmp = self.store.read_table(table_name)

        # Expand the data from unique
        if by_value and not (metadata.atom_metadata[property_name]["unique"]):
            tmp = self._build_atom_values(tmp, property_name)

        # Figure out unit scaling factors
        field_data = metadata.atom_metadata[property_name]
        if by_value and (field_data["units"] is not None) and (utype is not None):
            scale_factor = units.conversion_factor(field_data["utype"], utype)
            tmp[field_data["required_columns"]] *= scale_factor

        return tmp

    def add_atom_parameter(self, property_name, value, uid=None, utype=None):
        """
        Adds atom parameters to the Datalayer object

        Parameters
        ----------
        property_name : str
            The name of the atom property to be added
        value : float
            The value of the property to be added
        uid : int, optional
            The uid to assign to this parameterized term.
        utype : list of Pint units, options
            Custom units for this particular addition, otherwise uses the default units in the registered functional form.

        Results
        -------
        uid : int
            The set uid of the new atom parameter.

        Notes
        -----
        If a uid is not set and the parameter is already known, the function will return the current internal parameter.

        Examples
        --------

        assert 0 == dl.add_atom_parameter("charge", -0.8)
        assert 0 == dl.add_atom_parameter("charge", -0.8, uid=0)

        assert 5 == dl.add_atom_parameter("charge", -1.2, uid=5)
        assert 6 == dl.add_atom_parameter("charge", -2.e-19, uid=6, utype="coulomb")
        """

        property_name = self._check_atoms_dict(property_name)
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
            if found_key is not None:
                if found_key == uid:
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
        Lists all atom properties contained in the DataLayer.
        """
        return [k for k, v in self._atom_counts.items() if v > 0]

    def get_atom_count(self, property_name=None):
        """
        Get the number of atom properties added, or the maximum number of properties given.
        """

        # Get the current maximum
        if property_name is None:
            return max(v for k, v in self._atom_counts.items())

        property_name = self._check_atoms_dict(property_name)
        if property_name in self._atom_counts:
            return self._atom_counts[property_name]
        else:
            raise KeyError("DataLayer:get_atom_count: property_name `%s` not understood" % property_name)

    def list_atom_uids(self, property_name):

        property_name = self._check_atoms_dict(property_name)
        if not metadata.atom_metadata[property_name]["unique"]:
            return list(self._atom_metadata[property_name]["inv_uvals"])
        else:
            raise KeyError("DataLayere:list_atom_uids: '%s' is not stored as unique values." % property_name)

    def get_atom_parameter(self, property_name, uid, utype=None):
        """
        Obtains a atom parameter from the unique cache.
        """

        property_name = self._check_atoms_dict(property_name)
        if not metadata.atom_metadata[property_name]["unique"]:
            if not uid in self._atom_metadata[property_name]["inv_uvals"]:
                raise Exception("DataLayer:get_atom_parameter: property '%s' key '%d' not found." % (property_name,
                                                                                                     uid))
            cf = 1
            if utype is not None:
                cf = units.conversion_factor(metadata.atom_metadata[property_name]["utype"], utype)
            return self._atom_metadata[property_name]["inv_uvals"][uid] * cf
        else:
            raise KeyError("DataLayere:get_atom_parameter: '%s' is not stored as unique values." % property_name)

    def add_atoms(self, atom_df, by_value=False, utype=None):
        """
        Adds atom information to the DataLayer object.

        Parameters
        ----------
        atom_df : {DataFrame, list, tuple}
            The atom data to add to the object.
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

    def add_term_parameter(self, order, term_name, term_parameters, uid=None, utype=None):
        """
        Adds parameters for a given fuctional form.

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
        utype : list of Pint units, optional
            Custom units for this particular addition, otherwise uses the default units in the registered functional form.

        Return
        ------
        uid : ind
            The uid for the set functional form

        Examples
        --------

        assert 0 == dl.add_parameter(2, "harmonic", [4.0, 5.0])
        assert 0 == dl.add_parameter(2, "harmonic", [4.0, 5.0], uid=1)
        assert 1 == dl.add_parameter(2, "harmonic", [4.0, 6.0])

        """

        user_order = order
        order = metadata.sanitize_term_order_name(order)

        # Make sure we know what this is
        try:
            term_md = metadata.get_term_metadata(order, "forms", term_name)
        except KeyError:
            raise KeyError("DataLayer:add_parameters: Did not understand term order: %d, name: %s'." % (order,
                                                                                                        term_name))

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

    def get_term_parameter(self, order, uid, utype=None):

        order = metadata.sanitize_term_order_name(order)

        if uid not in self._terms[order]:
            raise KeyError("DataLayer:get_parameters: Did not find term '%d %d" % (order, uid))

        # Stored as [term_name, parameters...]
        data = self._terms[order][uid]

        term_md = metadata.get_term_metadata(order, "forms", data[0])

        # Zip up the parameters
        parameters = {k: v for k, v in zip(term_md["parameters"], data[1:])}

        # Were done
        if utype is None:
            return data[0], parameters

        # Need to convert
        if isinstance(utype, (list, tuple)):
            if len(utype) != len(term_md["parameters"]):
                raise KeyError("DataLayer:get_parameters: length of utype should match the length of parameters.")
            utype = {k: v for k, v in zip(term_md["parameters"], utype)}

        if not isinstance(utype, dict):
            raise TypeError("DataLayer:get_parameters: Input utype '%s' is not understood." % str(type(utype)))

        for key in term_md["parameters"]:
            parameters[key] *= units.conversion_factor(term_md["utype"][key], utype[key])

        return data[0], parameters

    def list_term_uids(self, order=None):

        # Return everything
        if order is None:
            ret = {}
            for k, v in self._terms.items():
                ret[k] = list(v)
            return ret

        # Just return a specific order
        order = metadata.sanitize_term_order_name(order)

        return list(self._terms[order])

    def get_term_count(self, order=None, uid=None):
        """
        Gives the number of term counts and uids

        Note
        ----
        The number of unique UID's may differ from list_parameter_uid's for incomplete DL's as
        `get_term_count` measures the number of terms add by `add_terms` while `list_parameter_uids` list
        the number of UID's added by `add_parameter`.
        """
        if (order is None) and (uid is None):
            return copy.deepcopy(self._term_count)

        if uid is None:
            if order is None:
                raise KeyError("DataLayer:get_term_count: Cannot use 'uid' if 'order' is None.")
            return self._term_count[order]

        return self._term_count[order][uid]

    def add_terms(self, order, df):
        """
        Adds terms using a index notation.

        Parameters
        ----------
        order : {str, int}
            The order (number of atoms) involved in the expression i.e. 2, "two"
        df : pd.DataFrame
            Adds a DataFrame containing the term information by index
            Required columns: ["term_index", "atom1_index", ..., "atom(order)_index", "term_index"]

        Returns
        -------
        return : bool
            Returns a boolean value if the operations was successful or not
        """

        order = metadata.sanitize_term_order_name(order)
        if order not in list(self._terms):
            raise KeyError("DataLayer:add_terms: Did not understand order key '%s'." % str(order))

        req_cols = metadata.get_term_metadata(order, "index_columns")

        not_found = set(req_cols) - set(df.columns)
        if not_found:
            raise KeyError("DataLayer:add_terms: Missing required columns '%s' for order %d" % (str(not_found), order))

        # Add the data in the index
        if "term_index" in df.columns:
            df = df[req_cols + ["term_index"]]
        else:
            raise Exception("NYI: Add terms by *not* term_index")

        # Get count information for later
        uvals, ucnts = np.unique(df["term_index"], return_counts=True)
        for uval, cnt in zip(uvals, ucnts):
            if uval not in self._term_count[order]:
                self._term_count[order][uval] = cnt
            else:
                self._term_count[order][uval] += cnt

            self._term_count[order]["total"] += cnt

        # Finally store the dataframe
        return self.store.add_table("term" + str(order), df)

    def get_terms(self, order):
        order = metadata.sanitize_term_order_name(order)
        if order not in list(self._terms):
            raise KeyError("DataLayer:add_terms: Did not understand order key '%s'." % str(order))

        try:
            return self.store.read_table("term" + str(order))
        except KeyError:
            cols = metadata.get_term_metadata(order, "index_columns") + ["term_index"]
            return pd.DataFrame(columns=cols)

    def add_bonds(self, bonds):
        """
        Adds bond using a index notation.

        Parameters
        ----------
        bonds : pd.DataFrame
            Adds a DataFrame containing the bond information by index
            Required columns: ["term_index", "atom1_index", "atom2_index", "term_type"]

        Returns
        -------
        return : bool
            Returns a boolean value if the operations was successful or not
        """

        return self.add_terms("bonds", bonds)

    def get_bonds(self):

        return self.get_terms("bonds")

    def add_angles(self, angles):
        """
        Adds angles using a index notation.

        Parameters
        ----------
        angles : pd.DataFrame
            Adds a DataFrame containing the angle information by index
            Required columns: ["term_index", "atom1_index", "atom2_index", "atom3_index", "term_type"]

        Returns
        -------
        return : bool
            Returns a boolean value if the operations was successful or not
        """

        return self.add_terms("angles", angles)

    def get_angles(self):

        return self.get_terms("angles")

    def add_dihedrals(self, dihedrals):
        """
        Adds dihedrals using a index notation.

        Parameters
        ----------
        dihedrals : pd.DataFrame
            Adds a DataFrame containing the dihedral information by index
            Required columns: ["term_index", "atom1_index", "atom2_index", "atom3_index", "atom4_index", "term_type"]

        Returns
        -------
        return : bool
            Returns a boolean value if the operations was successful or not
        """

        return self.add_terms("dihedrals", dihedrals)

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

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
from . import testing
from . import nb_converter

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

        # Create structure of _atom_metadata dictionary
        for k, v in metadata.atom_metadata.items():
            if not v["unique"]:
                self._atom_metadata[k] = {"uvals": {}, "inv_uvals": {}}
        self._atom_counts = {k: 0 for k in list(metadata.atom_metadata)}

        # Set up empty nonbond holders
        self._nb_parameters = {}
        self._nb_scaling_factors = {}
        self._nb_metadata = {}

        # Any remaining metadata
        self._box_size = {}
        self._box_center = {}
        self._mixing_rule = ''

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

    def set_mixing_rule(self, mixing_rule):
        """
        Store a mixing rule in the datalayer.

        Parameters:
        ------------------------------------
        mixing_rule: str
            Mixing rule to apply to calculate nonbonded parameters for pairs of atoms. Valid mixing rules are listed in
            nb_converter.LJ_mixing_functions 

        """
        if not isinstance(mixing_rule, str):
            raise TypeError("Validate mixing rule: %s is not a string" % mixing_rule)
        mixing_metadata = nb_converter.LJ_mixing_functions

        keys = mixing_metadata.keys()

        if mixing_rule not in keys:
            raise ValueError("Mixing rule type %s not found" % mixing_rule)

        self._mixing_rule = mixing_rule


    def get_mixing_rule(self):
        """ Retrieve the stored mixing rule from the datalayer. Returns a string """

        ret = copy.deepcopy(self._mixing_rule)
        return ret

    def set_box_center(self, box_center, utype=None):
        """
        Sets the center of the box.
        """

        # Get box metadata
        box_metadata = metadata.box_metadata
        dimensions = box_metadata["center"]

        # Make sure we have all keywords that define a simulation box
        for k in dimensions:
            if k.lower() not in box_center and k.upper() not in box_center:
                raise KeyError("Could not find key '%s'." % k)

        if utype is not None:
            if not isinstance(utype, dict):
                raise TypeError("Validate term dict: Unit type '%s' not understood" % str(type(utype)))

            # Convert to internal units
            for k, v in dimensions.items():
                internal = units.convert_contexts(v)
                cf = units.conversion_factor(utype[k], internal)
                self._box_center[k] = cf * box_center[k]

        else:
            for k, v in box_center.items():
                self._box_center[k] = v

    def get_box_center(self, utype=None):
        """
        Gets the overall size of the box for the datalayer
        """
        ret = copy.deepcopy(self._box_center)

        # Get information for internal representation of box dimensions
        box_metadata = metadata.box_metadata
        dimensions = box_metadata["center"]

        if utype is not None and ret:
            if not isinstance(utype, dict):
                raise TypeError("Validate term dict: Unit type '%s' not understood" % str(type(utype)))

            # Convert to internal units
            for k, v in dimensions.items():
                internal = units.convert_contexts(v)
                cf = units.conversion_factor(internal, utype[k])
                ret[k] *= cf

            return ret

        else:
            return ret

    def set_nb_scaling_factors(self, nb_scaling_factors):
        """
        Sets the exclusion information for the datalayer

        Parameters
        ------
        nb_scaling_factor: dict
            Format is:
        nb_scaling_factors = {
            "coul":{
                "scale12": "dimensionless",
                "scale13": "dimensionless",
                "scale14": "dimensionless",
            },
            "vdw":{
                "scale12": "dimensionless",
                "scale13": "dimensionless",
                "scale14": "dimensionless",
            }
        }
        """
        if not isinstance(nb_scaling_factors, dict):
            raise TypeError("Exclusion information cannot be validated as dictionary '%s'"% str(type(nb_scaling_factors)))

        exclusions_metadata = metadata.exclusions
        # Make sure we have all keywords 
        if nb_scaling_factors.keys() != exclusions_metadata.keys():
            raise KeyError("Not all exclusion keywords are imported")
    
        # Make sure scaling factors make sense
        for ok, ov in nb_scaling_factors.items():
            for k, v in ov.items():
                if v > 1.0 or v < 0.0:    
                    raise ValueError("Exclusion value outside bounds '%s'." % v)

        self._nb_scaling_factors = nb_scaling_factors

    def get_nb_scaling_factors(self):
        """
        Retrieves nonbonded scaling factors from datalayer.
        
        """
        ret = copy.deepcopy(self._nb_scaling_factors)

        return ret

    def set_nb_pair_interaction(self):
        """
        Set a special interaction between two particles
        :return:
        """
        return False

    def set_pair_scalings(self, scaling_df):
        """
        Set scaling factor for nonbond interaction between two atoms using multi-level indexing.

        Parameters:
        --------------------
            scaling_df: DataFrame
            Columns of the dataframe should be - 
                atom_index1: int
                    The atom index of the first atom
                atom_index2: int
                    The atom index of the second atom
                vdw_scale: float
                coul_scale: float

        Returns:
        -------------------
        Returns: bool
            True if successful
        """
        possible_columns = [y for x in metadata.additional_metadata.nb_scaling.values() for y in x]
        # Check the columns of the dataframe
        for col in scaling_df.columns:
            if col not in possible_columns:
                raise KeyError ("Column %s not recognized in set_pair_scalings." %(col))

        # Check to make sure atom_type1 and atom_type2 are set in dataframe
        for col in metadata.additional_metadata.nb_scaling["index"]:
            if col not in scaling_df.columns:
                raise KeyError("%s not found in scaling dataframe (set_pair_scalings)" %(col))

            if not np.issubdtype(scaling_df[col].dtype, int):
                raise TypeError("%s column is type %s. Should be integer" %(col, scaling_df[col].dtype) )
        
        # Make sure at least one scaling factor is set
        if len(scaling_df.columns) < 3:
            raise ValueError("No scaling factors set in set_pair_scalings")
        
        # Check that scalings are type float
        
        # Build multi-level indexer
        index = pd.MultiIndex.from_arrays([scaling_df["atom_index1"], scaling_df["atom_index2"]])

        for l in ["vdw_scale", "coul_scale"]:
            if l in scaling_df.columns:
                df = pd.Series(scaling_df[l].tolist(), index=index)
                self.store.add_table(l, df)

        return True

    def get_pair_scalings(self, nb_labels=["vdw_scale", "coul_scale"]):
        """
        Get scaling factor for nonbond interaction between two atoms
 
        Parameters
        ------------------------------------
            nb_labels: list
        
        Returns
        ------------------------------------
            pd.DataFrame
        """

        for k in nb_labels:
            if k not in metadata.additional_metadata.nb_scaling["data"]:
                raise KeyError("%s is not a valid nb_scale type" %(k))
        
        rlist = []
        rlabels = []

        for label in nb_labels:
            rlabels.append(label)
            rlist.append(self.store.read_table(label))
        
        ret = pd.concat(rlist, axis=1)
        ret.columns = rlabels

        return ret

    def build_scaling_list(self):
        """
        Build pair scalings based on parameters set in set_nb_scaling_factors.
        """
        scaling_factors = self.get_nb_scaling_factors()

        for k, v in scaling_factors.items():
            for scale, val in v.items():
                order = int(scale[-1])
                terms = self.get_terms(order)
                store_df = pd.DataFrame()

                store_df["atom_index1"] = terms["atom1"]
                store_df["atom_index2"] = terms["atom"+scale[-1]]
                store_df[k + "_scale"] = val
                
                if not store_df.empty and val is not 1.:
                    self.set_pair_scalings(store_df)

        return True

    def set_box_size(self, lattice_const, utype=None):
        """
        Sets the box lattice constants for the datalayer

        Inputs
        -----------------------------
        lattice_const: dict
            Dictionary containing box dimensions
            { 'a' : [length],
             'b' : [length],
             'c': [length],
             'alpha': [angle],
             'beta': [angle],
             'gamma': [angle],
             }

        """

        # Get box metadata
        box_metadata = metadata.box_metadata
        dimensions = box_metadata["dimensions"]

        # Make sure we have all keywords that define a simulation box
        for k in dimensions:
            if k.lower() not in lattice_const and k.upper() not in lattice_const:
                raise KeyError("Could not find key '%s'." % k)

        if utype is not None:
            if not isinstance(utype, dict):
                raise TypeError("Validate term dict: Unit type '%s' not understood" % str(type(utype)))

            # Convert to internal units
            for k, v in dimensions.items():
                internal = units.convert_contexts(v)
                cf = units.conversion_factor(utype[k], internal)
                self._box_size[k] = cf * lattice_const[k]

        else:
            for k, v in lattice_const.items():
                self._box_size[k] = v

    def get_box_size(self, utype=None):
        """
        Gets the overall size of the box for the datalayer
        """
        ret = copy.deepcopy(self._box_size)

        # Get information for internal representation of box dimensions
        box_metadata = metadata.box_metadata
        dimensions = box_metadata["dimensions"]

        if utype is not None and ret:
            if not isinstance(utype, dict):
                raise TypeError("Validate term dict: Unit type '%s' not understood" % str(type(utype)))

            # Convert to internal units
            for k, v in dimensions.items():
                internal = units.convert_contexts(v)
                cf = units.conversion_factor(internal, utype[k])
                ret[k] *= cf

            return ret

        else:
            return ret

    def evaluate(self, utype=None):
        """
        Evaluate the current state of the energy expression.
        """

        return energy_eval.evaluate_energy_expression(self, utype=utype)

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

            if isinstance(gb_idx, tuple):
                gb_dict = {k: v for k, v in zip(cols, gb_idx)}
            else:
                gb_dict = {cols[0]: gb_idx}

            gb_hash = utility.hash(gb_dict)

            # Update dictionary if necessary
            if gb_hash not in param_dict["uvals"]:

                # Bidirectional dictionary
                new_key = utility.find_lowest_hole(list(param_dict["inv_uvals"]))
                param_dict["uvals"][gb_hash] = new_key
                param_dict["inv_uvals"][new_key] = gb_dict

            # Grab the unique and set
            uidx = param_dict["uvals"][gb_hash]
            ret_df.loc[udf.index, property_name] = uidx

        return ret_df

    def _build_atom_values(self, df, property_name):
        """
        Expands the unique parameters using the built in property_name dictionary.
        """
        field_data = metadata.atom_metadata[property_name]
        param_dict = self._atom_metadata[property_name]

        cols = field_data["required_columns"]
        ret_df = pd.DataFrame(index=df.index, columns=cols, dtype=field_data["dtype"])
        for gb_idx, udf in df.groupby(cols):
            gb_dict = param_dict["inv_uvals"][gb_idx]
            for col in cols:
                ret_df.loc[udf.index, col] = gb_dict[col]

        return ret_df

    def _parse_atom_utype(self, property_name, utype):
        field_data = metadata.atom_metadata[property_name]
        if not isinstance(utype, dict):
            if property_name == "xyz":
                utype = {"X": utype, "Y": utype, "Z": utype}
            elif len(field_data["units"]) > 1:
                raise KeyError("Cannot infer utype without a dictionary.")
            else:
                utype = {field_data["required_columns"][0]: utype}
        return utype

    def _store_atom_table(self, table_name, df, property_name, by_value, utype):
        """
        Internal way to store atom tables
        """

        field_data = metadata.atom_metadata[property_name]

        # Figure out unit scaling factors
        if by_value and (field_data["units"] is not None) and (utype is not None):
            utype = self._parse_atom_utype(property_name, utype)
            cf = units.conversion_dict(utype, field_data["utype"])
            df = df[field_data["required_columns"]] * pd.Series(cf)

        # Make sure our ints are ints and not accidentally floats
        if np.issubdtype(field_data["dtype"], int):
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
            utype = self._parse_atom_utype(property_name, utype)
            cf = units.conversion_dict(field_data["utype"], utype)
            tmp[field_data["required_columns"]] *= pd.Series(cf)

        return tmp

    def add_atom_parameter(self, property_name, value, uid=None, utype=None, allow_duplicates=False):
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

        If a  uid is set and the parameter is already known, the function will
        return a KeyError unless allow_duplicates is set to True

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

        # Parse value
        if not isinstance(value, (int, float, dict)):
            raise TypeError("DataLayer:add_atom_parameter: Did not understand input type '%s'." % str(type(value)))

        if isinstance(value, (int, float)):
            req_cols = field_data["required_columns"]
            if len(req_cols) != 1:
                raise TypeError("DataLayer:add_atom_parameter: Expected %d values, only recieved one." % len(req_cols))
            value = {req_cols[0]: value}

        if (utype is not None) and not isinstance(utype, dict):
            utype = {field_data["required_columns"][0]: utype}

        # if (utype is not None) and (field_data["units"] is not None):
        tmp = {"parameters": field_data["required_columns"], "utype": field_data["utype"]}
        value = metadata.validate_term_dict(property_name, tmp, value, utype=utype)
        value = {k: v for k, v in zip(field_data["required_columns"], value)}

        # Round the floats
        if field_data["dtype"] == float:
            value = {k: round(v, field_data["tol"]) for k, v in value.items()}

        value_hash = utility.hash(value)

        # Check if we have this key
        found_key = None
        if value_hash in param_dict["uvals"]:
            found_key = param_dict["uvals"][value_hash]

        # Brand new uid, return next in sequence
        if uid is None:
            if found_key is not None:
                return found_key

            new_key = utility.find_lowest_hole(list(param_dict["inv_uvals"]))
            param_dict["uvals"][value_hash] = new_key
            param_dict["inv_uvals"][new_key] = value
            return new_key

        # We have a uid
        else:
            # Fine if it matches internally, otherwise throw
            if found_key is not None:
                if found_key == uid:
                    return found_key
                elif allow_duplicates is False:
                    raise KeyError(
                        "DataLayer:add_atom_parameters: Tried to add value %s, but found in uid (%d) and current keys (%d)"
                        % (value, uid, found_key))

            param_dict["inv_uvals"][uid] = value
            param_dict["uvals"][value_hash] = uid
            return uid

    def get_atom_parameter(self, property_name, uid, utype=None):
        """
        Obtains a atom parameter from the unique cache.
        """

        property_name = self._check_atoms_dict(property_name)
        if metadata.atom_metadata[property_name]["unique"]:
            raise KeyError("DataLayer:get_atom_parameter: '%s' is not stored as unique values." % property_name)

        if not uid in self._atom_metadata[property_name]["inv_uvals"]:
            raise KeyError("DataLayer:get_atom_parameter: property '%s' key '%d' not found." % (property_name, uid))
        field_data = metadata.atom_metadata[property_name]
        req_fields = field_data["required_columns"]

        data = copy.deepcopy(self._atom_metadata[property_name]["inv_uvals"][uid])

        # Handle utype
        if utype is not None:
            if field_data["units"] is None:
                raise TypeError(
                    "DataLayer:get_atom_parameter: property '%s' does not have units, but `utype` was passed in.")

            if not isinstance(utype, dict):
                if len(req_fields) > 1:
                    raise TypeError(
                        "DataLayer:get_atom_parameter: unit order can not be interpreted for more than one field.")
                utype = {req_fields[0]: utype}

            if set(field_data["utype"]) != set(utype):
                raise KeyError("DataLayer:get_atom_paramter: required units '%s' does not match input units '%s'." %
                               (str(list(field_data["utype"])), str(list(utype))))

            for k, v in field_data["utype"].items():
                data[k] *= units.conversion_factor(v, utype[k])

        if len(req_fields) == 1:
            return data[req_fields[0]]
        else:
            return data

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

    def get_bond_count(self):
        return len(self.get_bonds())

    def get_angle_count(self):
        return len(self.get_angles())

    def get_dihedral_count(self):
        return len(self.get_dihedrals())

    def get_unique_atom_types(self):
        return np.unique(self.get_atoms('atom_type'))

    def list_atom_uids(self, property_name):
        """
        Returns the unique values for a given property.
        """

        property_name = self._check_atoms_dict(property_name)
        if not metadata.atom_metadata[property_name]["unique"]:
            return list(self._atom_metadata[property_name]["inv_uvals"])
        else:
            raise KeyError("DataLayer:list_atom_uids: '%s' is not stored as unique values." % property_name)

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

    def list_valid_atom_properties(self):
        """
        Returns all possible atom properties which can be stored in the datalayer
        """
        return list(metadata.atom_property_to_column)

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
                    "DataLayer:add_term_parameter: uid keyword must be of type int, found type '%s'." % type(uid))

            # If we exist this could get dangerous
            if uid in self._terms[order]:
                old_param = self._terms[order][uid]
                match = (old_param[0] == term_name) and np.allclose(old_param[1:], params)
                if not match:
                    raise KeyError(
                        "DataLayer:add_term_parameter: uid already exists, but does not much current parameters.")
                else:
                    return uid

            else:
                params.insert(0, term_name)
                self._terms[order][uid] = params

                return uid

    def get_term_parameter(self, order, uid=None, utype=None):

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

    def list_term_parameters(self, order):
        """
        Gives information for all terms of specified order

        Parameters
        ----------
        order : int
            The order of the functional form (2, 3, 4, ...)

        Return
        ----------
        return : dict
            Returns dictionary of form
                { uid : [form_type, parameters] }
        """
        if order not in list(self._terms.keys()):
            raise KeyError("No terms with order %s exist" % order)

        return self._terms[order]

    def list_term_uids(self, order=None):
        """
        Lists all stored UID's in the datalayer.
        """

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

    def get_term_definition(self, order, uid):
        """
        Gives full definition for term of specified order and uid

        Parameters
        ----------
        order : int
            The order of the functional form (2, 3, 4, ...)
        uid: int
            The uid of

        Return
        ----------
        return : dict
            Returns dictionary of form
                { uid : [form_type, parameters] }
        """

        data = self._terms[order][uid]
        form = metadata.get_term_metadata(order, "forms", data[0])["form"]

        return (data[0], form)

    def summary(self):
        print("EEX DataLayer Object\n")

        print("System name: %s" % self.name)
        print("----------------------------------------------")

        # Print atom and topology information
        print("Atom Count:                 %d" % self.get_atom_count())
        print("Bond Count:                 %d" % self.get_bond_count())
        print("Angle Count:                %d" % self.get_angle_count())
        print("Dihedral Count:             %d" % self.get_dihedral_count())
        print("----------------------------------------------")

        # Print information about bond, angle, dihedral parameters
        print("Number of bond parameters:     %s" % len(self.list_term_uids()[2]))
        print("Number of angle parameters:    %s" % len(self.list_term_uids()[3]))
        print("Number of dihedral parameters: %s" % len(self.list_term_uids()[4]))

        print("----------------------------------------------")

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

### Non-bonded parameter

    def add_nb_parameter(self, atom_type, nb_name, nb_parameters, nb_model=None, atom_type2=None, utype=None):
        """
        Stores nb parameters in data layer as dictionary

        Stored dictionary has form:
            { (atom_type1, atom_type2):
                    'form' : nb_name
                    'parameters' :{
                            parameter_name1: nb_parameter,
                            parameter_name2: nb_parameter,
                            },
            }

        Parameters
        ----------
        atom_type : int
            The first atom type for the interaction
        nb_name: int
            The name of the functional form (eg - "LJ", "Buckingham")
        nb_parameters: dict
            Parameters for the functional form - should match nb_name and nb_form (if specified)
        nb_form: str
            (optional) - The form of the input (ex 'epsilon/sigma' with 'nb_name=LJ' indicates the input parameters are
            'epsilon' and 'sigma'
        atom_type2: int
            The atom type for the second interaction (optional)
        utype: dict
            Units of nb_parameters

        Returns
        -------------
        return: bool
              Returns True if successful
        """

        param_dict = {}
        param_dict['form'] = nb_name

        # Validate atom_type exists

        # Get functional form and ensure nb_parameters fit - maybe need to write function in validator.py
        try:
            form_md = metadata.get_nb_metadata(nb_name, model=nb_model)
        except KeyError:
            raise KeyError("DataLayer:add_parameters: Did not understand nonbond form: %d, name: %s'." % (nb_name,
                                                                                                          nb_form))
        parameters = form_md['parameters']

        # Validate input parameters against form
        if isinstance(nb_parameters, (list, tuple)):
            if len(parameters) == len(nb_parameters):
                param_dict['parameters'] = {k: v for k, v in zip(parameters, nb_parameters)}
            else:
                raise ValueError("Input number of parameters (%s) and number of form parameters (%s) do not match." %
                                 (len(nb_parameters), len(parameters)))
        elif isinstance(nb_parameters, (dict)):
            if set(nb_parameters.keys()) != set(parameters):
                raise ValueError("Incorrect parameters entered for nonbond form %s %s" % (nb_name, nb_form))
            else:
                param_dict['parameters'] = nb_parameters

        # Validate correct number of units are passed for parameters
        if utype is not None:
            if isinstance(utype, (list, tuple)):
                if len(utype) != len(form_md["utype"]):
                    raise ValueError("Validate term dict: Number of units passed is %d, expected %d" %
                                     (len(utype), len(form["utype"])))
                form_units = list(utype)
            elif isinstance(utype, dict):
                form_units = []
                for key in parameters:
                    try:
                        form_units.append(utype[key])
                    except KeyError:
                        raise KeyError("Validate term dict: Did not find expected key '%s' from term'." % (key))
            else:
                raise TypeError("Validate term dict: Unit type '%s' not understood" % str(type(utype)))

            # Convert to internal units
            for x, key in enumerate(form_md["parameters"]):
                cf = units.conversion_factor(form_units[x], form_md["utype"][key])
                param_dict['parameters'][key] *= cf

        if (nb_name == "LJ"):
            model_default = metadata.get_nb_metadata(nb_name, "default")
            param_dict['parameters'] = nb_converter.convert_LJ_coeffs(param_dict['parameters'], nb_model,
                                                                      model_default)

        # Store it! --
        param_dict_key = (atom_type, atom_type2)

        # Sort if atom_type2 is not None. atom_type1 is always less than atom_type2
        if atom_type2 != None:
            param_dict_key = tuple(sorted(param_dict_key))


        self._nb_parameters[param_dict_key] = param_dict
        return True

    def get_nb_parameter(self, atom_type, nb_model=None, atom_type2=None, utype=None):
        """
        Retrieves nb parameter from datalayer

        Parameters
        -----------------
        atom_type: int
            The first atom_type of the NB interaction
        nb_form: str
            The desired output form (optional). If not indicated, default for datalayer will be returned.
        atom_type2: int
            The second atom_type of the NB interaction (optional). If atom_type and atom_type2 are specified, returned
            parameters apply to the interaction between atom_type1 and atom_type2
        utype: dict
            Units of output. Must be compatible with parameters.

        Returns
        ------------------
        return: dict

        Returned dictionary has form:
            {
                parameter_name_1: nb_parameter_1,
                parameter_name_2: nb_parameter_2,
                ...
                parameter_name_n : nb_parameter_n,
            }
        """

        # Build key
        # Don't need if statement here, just use (atom_type, None)
        param_dict_key = (atom_type, atom_type2)

        # Get information from data layer - check that interaction is set for atom types
        if param_dict_key in self._nb_parameters.keys():
            # Use deep copy here
            nb_parameters = copy.deepcopy(self._nb_parameters[param_dict_key])
        else:
            raise KeyError("Nonbond interaction for atom types (%s, %s) not found" % param_dict_key)

        # Get nb_name - this is stored when parameter is input (ex "LJ")
        nb_name = nb_parameters['form']

        # Grab data we want from data layer
        param_dict = nb_parameters["parameters"]

        # Get and validate datalayer units for nb form parameters (form is from metadata)
        form_md = metadata.get_nb_metadata(nb_name, model=nb_model)

        # Find models
        default_form = metadata.get_nb_metadata(nb_parameters["form"], "default")
        if nb_model is None:
            nb_model = default_form

        ### Need to convert to specified nb_name (form) if needed (ex - AB to epsilon/sigma)
        if nb_parameters["form"] == "LJ":
            param_dict = nb_converter.convert_LJ_coeffs(param_dict, metadata.get_nb_metadata("LJ", "default"),
                                                        nb_model)

        # Convert units if specified - otherwise return what is stored in datalayer
        form_units = {}
        if utype is not None:
            for key in form_md["parameters"]:
                try:
                    form_units[key] = utype[key]
                except KeyError:
                    raise KeyError("Validate term dict: Did not find expected key '%s' from term (utype)'." % (key))

            for x, key in enumerate(param_dict):
                # Convert from what is in DL (form["utype"][key] to user specified units (form_units[x]
                cf = units.conversion_factor(form_md["utype"][key], form_units[key])
                param_dict[key] *= cf

        return param_dict

    def mix_LJ_parameters(self, atom_type1, atom_type2, mixing_rule=None):
        """
        Mixes LJ parameters based on atom types and mixing rules. Stores in datalayer as (atom_type1, atom_type2)

        Parameters:
        ---------------------------------
        atom_type: int
            The first atom_type of the NB interaction
        atom_type: int
            The second atom type of the NB interaction
        mixing_rule: str
            The mixing rule to combine parameters with. If not set, will be set to mixing rule stored in datalayer


        Returns:
        --------------------------------
        return: bool
            Returns True if successful
        """

        if mixing_rule == None:
            mixing_rule = self._mixing_rule

        param_keys = [(atom_type1, None), (atom_type2, None)]

        params = []

        # Get information from data layer - check that interaction is set for atom types
        for k in param_keys:
            if k in self._nb_parameters.keys():
                # Use deep copy here
                params.append(copy.deepcopy(self._nb_parameters[k]))
            else:
                raise KeyError("Nonbond interaction for atom types (%s, %s) not found" % (k))

        # Check that both parameters are LJ form
        if params[0]["form"] != "LJ" or params[1]["form"] != "LJ":
            raise ValueError("Can only combine LJ coefficients using mixing rules.")

        # Apply mixing rule
        new_params = nb_converter.mix_LJ(params[0]["parameters"], params[1]["parameters"], mixing_rule=mixing_rule)

        # Add new parameter to datalayer!
        self.add_nb_parameter(atom_type=atom_type1, atom_type2=atom_type2, nb_parameters=new_params,
                              nb_name="LJ", nb_model="AB")

        return True

    def build_LJ_mixing_table(self):
        """
        Function applies mixing rule to (atom_type, None) pairs.

        Stores NB parameters in datalayer dictionary with keys
            (1...n, 1...n)
             where n is the number of atom types.

        Returns: bool
            Returns True if successful
        """

        parameters = self.list_nb_parameters(nb_name="LJ", nb_model="AB", itype="single")

        parameter_keys = list(parameters)

        for k in range(0, len(parameter_keys)):
            for k2 in range(0, k+1):
                self.mix_LJ_parameters(atom_type1=parameter_keys[k][0], atom_type2=parameter_keys[k2][0])


        return True


    def list_stored_nb_types(self):
        """
        Returns the type of NB interactions stored in datalayer (ex Lennard Jones, or Buckingham) as a list
            ex. - ["LJ", "Buckingham"]
        """

        nb_types = []
        for k, v in self._nb_parameters.items():
            nb_types.append(v['form'])
        unique_nb_types = np.unique(nb_types)
        return unique_nb_types

    def list_nb_parameters(self, nb_name, nb_model=None, utype=None, itype="all"):
        """
        Return all NB parameters stored in data layer which have the form specified by nb_name.

        Parameters
        ------------------
        nb_name: str
            Name of nonbond potential (ex "LJ" or "Buckingham") to list interactions for
        nb_form: str (optional)
            Output form of potential (ex "epsilon/sigma" for nb_name "LJ"). If not specified, default for datalayer will
            be returned
        utype: dict (optional)
            Units for output. Must be compatible with nb_name and form. If not specified, default for datalayer will be
            returned
        itype: interaction type (optional)
            Can specify "all", "pair", or "single". All returns all values in the datalayer, while pair returns interactions
            for I J pairs (i.e. atom_type1, atom_type2) and "single" returns only (atom_type1, None) interactions.

        Returns
        ------------------
        return_parameters: dict
            Returned dict has form
                { (atom_type1, atom_type_2)
                        {
                            'nb_parameter_name' : value,
                            'nb_parameter_name2' : value,
                        }
                }
        """
        return_parameters = {}
        term_dict = self._nb_parameters.copy()

        for key, value in term_dict.items():
            if value['form'] != nb_name:
                continue

            if (None in key) and (itype == "all" or itype == "single"):
                return_parameters[key] = self.get_nb_parameter(
                    atom_type=key[0], atom_type2=key[1], nb_model=nb_model, utype=utype)
            elif (None not in key) and (itype == "all" or itype == "pair"):
                return_parameters[key] = self.get_nb_parameter(
                    atom_type=key[0], atom_type2=key[1], nb_model=nb_model, utype=utype)

        return return_parameters




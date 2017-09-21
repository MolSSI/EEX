"""
Contains the DataLayer class (name in progress) which takes and reads various peices of data
"""

import os
import pandas as pd
import numpy as np
import tables
import collections

from . import filelayer
from . import metadata
from . import units


APC_DICT = metadata.atom_property_to_column

class DataLayer(object):
    def __init__(self, name, store_location=None, save_data=False, backend="HDF5"):
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
        self.parameters = {}
        self.terms = {2: {}, 3: {}, 4:{}, 5:{}, 6:{}, 7:{}}


    def _validate_table_input(self, data, needed_cols):
        if isinstance(data, (list, tuple)):
            if len(data) != len(needed_cols):
                raise Exception("DataLayer:add_data requires exactly the %s columns" % str(needed_cols))
            data = pd.DataFrame([data], columns=needed_cols)
        elif isinstance(data, pd.DataFrame):
            if collections.Counter(needed_cols) != collections.Counter(data.columns):
                raise Exception("DataLayer:add_data requires exactly the %s columns" % str(needed_cols))
        else:
            raise TypeError("DataLayer:add_data type %s not recognized" % type(data))

        return data

    def _set_unique_params(self, df, parameter_name):
        """
        Hashes the input parameters to build in internal index of unique values.
        """

        field_data = metadata.atom_metadata[parameter_name]
        if parameter_name not in list(self.parameters):
            self.parameters[parameter_name] = {"counter": -1, "uvals": {}, "inv_uvals": {}}

        cols = field_data["required_columns"]

        if field_data["dtype"] == float:
            df = df[cols].round(field_data["tol"])

        ret_df = pd.DataFrame(index=df.index)
        ret_df[parameter_name] = 0
        param_dict = self.parameters[parameter_name]


        # For each unique value
        for gb_idx, udf in df.groupby(cols):

            # Update dictionary if necessary
            if not gb_idx in list(param_dict["uvals"]):
                param_dict["counter"] += 1

                # Bidirectional dictionary
                param_dict["uvals"][gb_idx] = param_dict["counter"]
                param_dict["inv_uvals"][param_dict["counter"]] = gb_idx

            # Grab the unique and set
            uidx = param_dict["uvals"][gb_idx]
            ret_df.loc[udf.index, parameter_name] = uidx

        return ret_df

    def _build_value_params(self, df, parameter_name):
        """
        Expands the unique parameters using the built in parameter_name dictionary.
        """
        field_data = metadata.atom_metadata[parameter_name]
        param_dict = self.parameters[parameter_name]

        cols = field_data["required_columns"]

        ret_df = pd.DataFrame(index=df.index)
        ret_df[parameter_name] = 0.0

        for gb_idx, udf in df.groupby(cols):
            ret_df.loc[udf.index, parameter_name] = param_dict["inv_uvals"][gb_idx]

        return ret_df

    def _store_table(self, table_name, df, parameter_name, by_value):

        if by_value and not (metadata.atom_metadata[parameter_name]["unique"]):
            tmp_df = self._set_unique_params(df, parameter_name)
        else:
            tmp_df = df[APC_DICT[parameter_name]]

        return self.store.add_table(table_name, tmp_df)

    def _get_table(self, table_name, parameter_name, by_value):

        tmp = self.store.read_table(table_name)
        if by_value and not (metadata.atom_metadata[parameter_name]["unique"]):
            tmp = self._build_value_params(tmp, parameter_name)
        return tmp


    def add_atoms(self, atom_df, property_name=None, by_value=False):
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

        Returns
        -------
        return : {bool}
            If the add was successful or not.


        Example
        -------
        dl = DataLayer("test")

        # Add the atoms to the same molecule
        dl.add_atoms([0, 1], properties="molecule_index")
        dl.add_atoms([1, 1], properties="molecule_index")
        dl.add_atoms([2, 1], properties="molecule_index")

        # Add the XYZ information for five random atoms
        tmp_df = pd.DataFrame(np.random.rand(5, 3), columns=["X", "Y", "Z"])
        tmp_df["atom_index"] = np.arange(5)
        dl.add_atom(tmp_df)
        """

        # Our index name
        index = "atom_index"
        if property_name and (property_name not in list(APC_DICT)):
            raise KeyError("DataLayer:add_atoms: Property name '%s' not recognized." % property_name)

        # Parse list and DataFrame object and check validity
        if isinstance(atom_df, (list, tuple)):
            if property_name is None:
                raise KeyError("DataLayer:add_atoms: If data type is list, 'property_name' must be set.")
            if len(atom_df) != (len(APC_DICT[property_name]) + 1):
                raise KeyError("DataLayer:add_atoms: Property name '%s' not recognized." % property_name)

            atom_df = pd.DataFrame([atom_df], columns=[index] + APC_DICT[property_name])
            atom_df.set_index(index, drop=True, inplace=True)

        elif isinstance(atom_df, pd.DataFrame):
            if index in atom_df.columns:
                atom_df = atom_df.set_index(index, drop=True)

            if atom_df.index.name != index:
                raise KeyError("DataLayer:add_atoms: DF index must be the `atom_index`.")
        else:
            raise KeyError("DataLayer:add_atoms: Data type '%s' not understood." % type(atoms_df))

        # Add a single property
        if property_name:
            self._store_table(property_name, atom_df, property_name, by_value)
        # Try to add all possible properties
        else:
            set_cols = set(atom_df.columns)
            found_one = False
            for k, v in APC_DICT.items():
                # Check if v is in the set_cols (set logic)
                if set(v) <= set_cols:
                    self._store_table(k, atom_df, k, by_value)
                    found_one = True
            if not found_one:
                raise Exception("DataLayer:add_atom: No data was added as no key was matched from input columns:\n%s" %
                                (" " * 11 + str(atom_df.columns)))

        return True

    def get_atoms(self, properties, by_value=False):
        """
        Obtains atom information to the DataLayer object.

        Parameters
        ----------
        properties : {list, str}
            The properties to obtain for the atom data.

        Returns
        -------
        return : pd.DataFrame
            Returns a DataFrame containing the atom property information
            If the add was successful or not.

        """

        valid_properties = list(metadata.atom_property_to_column)

        # Our index name
        index = "atom_index"
        if not isinstance(properties, (tuple, list)):
            properties = [properties]

        if not set(properties) <= set(list(valid_properties)):
            invalid_props = set(properties) - set(list(valid_properties))
            raise KeyError("DataLayer:add_atoms: Property name(s) '%s' not recognized." % str(list(invalid_props)))

        df_data = []
        for prop in properties:
            # tmp = self.store.read_table(prop)
            # if by_value:
            #     print(tmp.head())
            #     print(prop)
            #     tmp = self._build_value_params(tmp, prop)
            tmp = self._get_table(prop, prop, by_value)
            df_data.append(tmp)
            # df_data.append(self.store.read_table(prop))

        return pd.concat(df_data, axis=1)

    def add_functional_form(self, order, functional_form, term_parameters, units=None):
        """
        Adds a n-body energy expression to the DataLayer object.
        """

        term_md = metadata.get_term_data(order, functional_form)

        if len(term_parameters) != (term_md["parameters"]):
            raise ValueError("DataLayer:add_functional_form: Length of parameters is %d, expected %d for term '%s'." %
                             (len(term_parameters), len(term_md["parameters"]), functional_form))

        # Figure out if this a new term or not
        if units:
            if len(term_parameters) != len(units):
                raise ValueError("DataLayer:add_functional_form: Length of input parameters and units must match.")
            params = [units.ureg.Quantity(v, u) for v, u in zip(term_parameters, units)]

        self.terms[order]

        # index_1, index_2, term_uuid
        # DL.register_functional_forms(metadata)
        # uuid = DL.add_parameters(2, "class2", [0, 3, 4, 2], uuid=xx)
        # DL.add_term(2, [[index_1, index_2, uuid]])
        # DL.add_bond([[index_1, index_2, uuid]])

        # DL.add_term([[index_1, index_2, class, [0, 3, 4, 2]])V


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

        needed_cols = ["bond_index", "atom1_index", "atom2_index", "bond_type"]

        bonds = self._validate_table_input(bonds, needed_cols)

        # Reorder columns
        bonds = bonds[needed_cols]

        self.store.add_table("bonds", bonds)

        return True

    def get_bonds(self):

        return self.store.read_table("bonds")

    def add_angles(self, angles):
        """
        Adds angles using a index notation.

        Parameters
        ----------
        bonds : pd.DataFrame
            Adds a DataFrame containing the angle information by index
            Required columns: ["bond_index", "atom1_index", "atom2_index", "atom3_index", "bond_type"]

        Returns
        -------
        return : bool
            Returns a boolean value if the operations was successful or not
        """

        needed_cols = ["angle_index", "atom1_index", "atom2_index", "atom3_index", "angle_type"]

        angles = self._validate_table_input(angles, needed_cols)

        # Reorder columns
        angles = angles[needed_cols]

        self.store.add_table("angles", angles)

    def get_angles(self):

        return self.store.read_table("angles")

    def call_by_string(self, *args, **kwargs):

        if args[0] == "NYI":
            return

        try:
            function = getattr(self, args[0])
        except:
            raise KeyError("DataLayer:call_by_string: does not have method %s." % args[0])

        return function(*args[1:], **kwargs)

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

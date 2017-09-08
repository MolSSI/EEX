"""
Contains the DataLayer class (name in progress) which takes and reads various peices of data
"""

import os
import pandas as pd
import tables
import collections

from . import filelayer
from . import fields


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

    def add_atoms(self, atom_df, property_name=None):
        """
        Adds atom information to the DataLayer object.

        Parameters
        ----------
        atom_df : {DataFrame, list, tuple}
            The atom data to add to the object.
        property_name: {list, str}, optional
            The atom property that is added, only necessary if a list is passed in.

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

        apc_dict = fields.atom_property_to_column

        # Our index name
        index = "atom_index"
        if property_name and (property_name not in list(apc_dict)):
            raise KeyError("DataLayer:add_atoms: Property name '%s' not recognized." % property_name)

        # Parse list and DataFrame object and check validity
        if isinstance(atom_df, (list, tuple)):
            if property_name is None:
                raise KeyError("DataLayer:add_atoms: If data type is list, 'property_name' must be set.")
            if len(atom_df) != (len(apc_dict[property_name]) + 1):
                raise KeyError("DataLayer:add_atoms: Property name '%s' not recognized." % property_name)

            atom_df = pd.DataFrame([atom_df], columns=[index] + apc_dict[property_name])
            atom_df.set_index(index, drop=True, inplace=True)

        elif isinstance(atom_df, pd.DataFrame):
            if index in atom_df.columns:
                atom_df = atom_df.set_index(index, drop=True)
            else:
                atom_df.index.name = index
        else:
            raise KeyError("DataLayer:add_atoms: Data type '%s' not understood." % type(atoms_df))

        # Add a single property
        if property_name:
            self.store.add_table(property_name, atom_df[apc_dict[property_name]])
        # Try to add all possible properties
        else:
            set_cols = set(atom_df.columns)
            found_one = False
            for k, v in apc_dict.items():
                # Check if v is in the set_cols (set logic)
                if set(v) <= set_cols:
                    self.store.add_table(k, atom_df[v])
                    found_one = True
            if not found_one:
                raise Exception(
                    "DataLayer:add_atom: No data was added as no key was matched from input columns:\n%s" %
                    (" " * 11 + str(atom_df.columns)))

        return True

    def get_atoms(self, properties):
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

        valid_properties = list(fields.atom_property_to_column)

        # Our index name
        index = "atom_index"
        if not isinstance(properties, (tuple, list)):
            properties = [properties]

        if not set(properties) <= set(list(valid_properties)):
            invalid_props = set(properties) - set(list(valid_properties))
            raise KeyError("DataLayer:add_atoms: Property name(s) '%s' not recognized." % str(list(invalid_props)))

        df_data = []
        for prop in properties:
            df_data.append(self.store.read_table(prop))

        return pd.concat(df_data, axis=1)

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

        needed_cols = ["bond_index", "atom1_index", "atom2_index", "atom3_index", "angle_type"]

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

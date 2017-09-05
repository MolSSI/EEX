"""
Contains the DataLayer class (name in progress) which takes and reads various peices of data
"""

import os
import pandas as pd
import tables
import collections

from . import filelayer


# Valid Atom Properties
_valid_atom_properties = {
    "molecule_index": ["molecule_index"],
    "atom_type": ["atom_type"],
    "charge": ["charge"],
    "XYZ": ["X", "Y", "Z"]
}

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
        """

        # Set the state
        self.name = name

        # Build the store
        self.store_location = store_location
        if self.store_location is None:
            self.store_location = os.getcwd()

        if backend == "HDF5":
            self.store = filelayer.HDFStore(self.name, self.store_location, save_data)
        else:
            raise KeyError("DataLayer:add_atoms ")

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

        # Our index name
        index = "atom_index"
        if property_name and (property_name not in list(_valid_atom_properties)):
            raise KeyError("DataLayer:add_atoms: Property name '%s' not recognized." % property_name)

        if isinstance(atom_df, (list, tuple)):
            if property_name is None:
                raise KeyError("DataLayer:add_atoms: If data type is list, 'property_name' must be set.")
            if len(atom_df) != (len(_valid_atom_properties[property_name]) + 1):
                raise KeyError("DataLayer:add_atoms: Property name '%s' not recognized." % property_name)

            atom_df = pd.DataFrame([atom_df], columns=[index] + _valid_atom_properties[property_name])
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
            self.store.add_table(property_name, atom_df[_valid_atom_properties[property_name]])
        # Try to add all possible properties
        else:
            set_cols = set(atom_df.columns)
            for k, v in _valid_atom_properties.items():
                if set(v) <= set_cols:
                    self.store.add_table(k, atom_df[v])

    def get_atoms(self, properties):
        # Our index name
        index = "atom_index"
        if not isinstance(properties, (tuple, list)):
            properties = [properties]

        if not set(properties) <= set(list(_valid_atom_properties)):
            invalid_props = set(properties) - set(list(_valid_atom_properties))
            raise KeyError("DataLayer:add_atoms: Property name(s) '%s' not recognized." % str(list(invalid_props)))

        df_data = []
        for prop in properties:
            df_data.append(self.store.read_table(prop))

        return pd.concat(df_data, axis=1)

    def add_bonds(self, bonds):
        """
        Adds bond using a index notation:

        Bonds:

        Parameters
        ----------
        """

        needed_cols = ["bond_index", "atom1_index", "atom2_index", "bond_type"]

        bonds = self._validate_table_input(bonds, needed_cols)

        # Reorder columns
        bonds = bonds[needed_cols]

        self.store.add_table("bonds", bonds)

    def get_bonds(self):

        return self.store.read_table("bonds")

    def add_angles(self, angles):
        """
        Adds bond using a index notation:

        Bonds:

        Parameters
        ----------
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

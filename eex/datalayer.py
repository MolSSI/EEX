"""
Contains the DataLayer class (name in progress) which takes and reads various peices of data
"""

import os
import pandas as pd
import tables
import collections

from . import filelayer


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

    def add_atoms(self, atoms):

        needed_cols = ["atom_index", "molecule_index", "atom_type", "X", "Y", "Z"]

        if isinstance(atoms, (list, tuple)):
            if len(atoms) != len(needed_cols):
                raise Exception("DataLayer:add_atoms requires exactly the %s columns" % str(needed_cols))
            atoms = pd.DataFrame([atoms], columns=needed_cols)
        elif isinstance(atoms, pd.DataFrame):
            if collections.Counter(needed_cols) != collections.Counter(atoms.columns):
                raise Exception("DataLayer:add_atoms requires exactly the %s columns" % str(needed_cols))
        else:
            raise TypeError("DataLayer:add_atoms type %s not recognized" % type(atoms))

        # Reorder columns
        atoms = atoms[needed_cols]

        # Store the Data
        self.store.add_table("atoms", atoms)

    def read_atoms(self):

        return self.store.read_table("atoms")

    def add_bonds_by_index(self, bonds):
        """
        Adds bond using a index notation:

        Bonds:




        Parameters
        ----------
        """
        raise Exception("NYI")

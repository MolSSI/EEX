"""
Base class for the filelayer.
"""

import os
import pandas as pd
import tables

class HDFStore(object):
    def __init__(self, name, store_location, save_data):


        # Figure out what the store_location means
        self.name = name
        if store_location is None:
            self.store_location = os.getcwd()
        else:
            self.store_location = store_location

        if not os.path.exists(self.store_location):
            os.makedirs(self.store_location)

        # Setup the store
        self.save_data = save_data
        self.store_filename = os.path.join(self.store_location, self.name + ".h5")
        self.store = pd.HDFStore(self.store_filename)

        # Set additional state
        self.created_tables = []


    def add_table(self, key, data):
        """
        Appends or builds a new table

        Parameters
        ----------
        key : str
            The name of the table
        data : pd.DataFrame
            The data to append to the table
        """

        # Do we append, or do we need a new table?
        do_append = True
        if key not in self.created_tables:
            do_append = False
            self.created_tables.append(key)

        data.to_hdf(self.store, key, format="t", append=do_append)

    def read_table(self, key, rows=None):

        if rows:
            raise Exception("NYI")

        return pd.read_hdf(self.store, key)

    def __del__(self):
        """
        On objection deletion close the store and remove store (optional)
        """

        self.store.close()
        if not self.save_data:
            os.unlink(self.store_filename)

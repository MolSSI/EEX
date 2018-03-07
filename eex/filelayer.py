"""
Base class for the filelayer.
"""

import os

import pandas as pd


def build_store(store_type, name, store_location, save_data):

    if store_type.upper() == "HDF5":
        return HDFStore(name, store_location, save_data)
    elif store_type.upper() == "MEMORY":
        return MemoryStore(name, store_location, save_data)
    else:
        raise KeyError("build_store: store_type of type '%s' not recognized." % store_type)


class BaseStore(object):
    def __init__(self, name, store_location, save_data):

        # Figure out what the store_location means
        self.name = name
        self.store_location = self._parse_store_location(store_location)

        # Save other flags
        self.save_data = save_data

    def _parse_store_location(self, store_location):
        if store_location is None:
            store_location = os.getcwd()
        else:
            store_location = store_location

        if not os.path.exists(store_location):
            os.makedirs(store_location)

        return store_location


class HDFStore(BaseStore):
    def __init__(self, name, store_location, save_data):

        # Init the base class
        BaseStore.__init__(self, name, store_location, save_data)

        # Setup the store
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

        # We do nothing for zero data
        if 0 in data.shape:
            return False

        # Do we append, or do we need a new table?
        do_append = True
        if key not in self.created_tables:
            do_append = False
            self.created_tables.append(key)

        data.to_hdf(self.store, key, format="t", append=do_append)
        return True

    def read_table(self, key, rows=None, where=None, chunksize=None):
        """
        Reads the table using either the rows or where syntax
        Appends or builds a new table

        Parameters
        ----------
        key : str
            The name of the table
        data : pd.DataFrame
            The data to append to the table
        """

        if rows:
            raise Exception("NYI")
        
        if key in self.list_tables():
            return pd.read_hdf(self.store, key)
        else:
            return pd.DataFrame()

    def close(self):
        """
        Closes the FL file.
        """

        self.store.close()
        if not self.save_data:
            try:
                os.unlink(self.store_filename)
            except OSError:
                pass

    def list_tables(self):

        return list(self.created_tables)

    def __del__(self):
        """
        On objection deletion close the store and remove store (optional)
        """

        self.close()

    def copy_table(self, from_key, to_key, columns_rename=None):
        """
        Copies a table from one key to another
        """

        for chunk in pd.read_hdf(self.store, from_key, iterator=True, chunksize=1.e6):
            if columns_rename is not None:
                chunk.rename(columns=columns_rename, inplace=True)
            self.add_table(to_key, chunk)


class MemoryStore(BaseStore):
    def __init__(self, name, store_location, save_data):

        # Init the base class
        BaseStore.__init__(self, name, store_location, save_data)
        self.store_filename = os.path.join(self.store_location, self.name + ".h5")

        # Table holder dictionary
        self.tables = {}
        self.table_frags = {}

    def add_table(self, key, data):

        # Lazy concat fragments for speed
        if key not in list(self.tables):
            self.table_frags[key] = [data]
            self.tables[key] = pd.DataFrame()
        else:
            self.table_frags[key].append(data)

    def read_table(self, key):
        if key not in list(self.tables):
            raise KeyError("Key %s does not exist" % key)

        # Concat any fragments
        if len(self.table_frags[key]):
            if not self.tables[key].empty:
                self.table_frags[key].insert(0, self.tables[key]) 
            self.tables[key] = pd.concat(self.table_frags[key])
            self.table_frags[key] = []

        return self.tables[key].copy()

    def close(self):
        """
        Closes the FL file.
        """

        if self.save_data:
            store = pd.HDFStore(self.store_filename)
            for k, v in self.tables.items():
                v.to_hdf(store, k, format="t")
            store.close()

    def list_tables(self):

        return list(self.tables)

    def copy_table(self, from_key, to_key, columns_rename=None):
        """
        Copies a table from one key to another
        """

        tmp = self.read_table(from_key)

        if columns_rename is not None:
            tmp.rename(columns=columns_rename, inplace=True)

        self.add_table(to_key, tmp)

    def __del__(self):
        self.close()

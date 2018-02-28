# EEX Structure

## Overview of EEX Data Storage Strategy

Data in EEX is stored in two different places: the `FileLayer` and the `DataLayer`.  

The `FileLayer` data can either be stored in memory or on disk as HDF5. Thus, the `FileLayer` is for storing information which has the potential to be large and heterogenous. Currently, the file layer contains data for each atom, and topology information for the system (bonds, angles, dihedrals). 

The `DataLayer` is for storing smaller, less uniform data. Force field parameters are stored in the datalayer 

## Private Data Structures


### Atoms

Information for each individual atom is stored in the file layer in store.??

#### dl._atom_metadata
`dict`
This dictionary is used for atom properties which are *not* unique (i.e. several atoms will have the same values).    
Has following structure:  

``` 
atom_metadata = {
		'atom_property1' :
			{ 
				'uvals' : {
					    uid : value, 
					  }
				'inv_uvals' : {
					    value: uid
						},
		'atom_property2' : ...
		
		}
```

where `atom_property1`, `atom_property2` come from atom_metadata.py  

This is used in functions like `add_atom_parameter`, `get_atom_parameter`, `list_atom_uids` (user accessible), and in internal functions `_find_unique_atom_values`,`_build_atom_values`, which map properties stored in this way (through `add_atom_parameter`)
rather than `add_atoms`, for example, to stored atom information. When this is done, it is assumed that they uid corresponds to the atom_type.

How does atom metadata work in the datalayer?
add_atoms
_store_atom_tabl


#### dl._atom_counts
dict
Lists the number of values stored for each metadata item

### Bonded Terms

#### dl._terms
dict

#### dl._term_count
dict

### Nonbonded Terms



# File Layer

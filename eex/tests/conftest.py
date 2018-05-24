import eex
import numpy as np
import pytest
import pandas as pd
import copy
from . import eex_find_files

@pytest.fixture(scope="function", params=["HDF5", "Memory"])
def butane_dl(request):
    # Build the topology for UA butane. Force Field parameters can be set to defaults by using
    # ff=True, nb=True, or scale=True (defaults)
    dl = eex.datalayer.DataLayer("butane", backend=request.param)

    def build_butane(ff=True, nb=True, scale=True):

        # Create empty dataframe
        atom_df = pd.DataFrame()

        # Create atomic system using pandas dataframe
        atom_df["atom_index"] = np.arange(0, 4)
        atom_df["molecule_index"] = [int(x) for x in np.zeros(4)]
        atom_df["residue_index"] = [int(x) for x in np.zeros(4)]
        atom_df["atom_name"] = ["C1", "C2", "C3", "C4"]
        atom_df["charge"] = np.zeros(4)
        atom_df["atom_type"] = [1, 2, 2, 1]
        atom_df["X"] = [0, 0, 0, -1.474]
        atom_df["Y"] = [-0.4597, 0, 1.598, 1.573]
        atom_df["Z"] = [-1.5302, 0, 0, -0.6167]
        atom_df["mass"] = [15.0452, 14.02658, 14.02658, 15.0452]

        # Add atoms to datalayer
        dl.add_atoms(atom_df, by_value=True)

        # Create empty dataframes for bonds
        bond_df = pd.DataFrame()

        # Create column names. Here, "term_index" refers to the bond type index.
        # i.e. - if all bonds are the same type, they will have the same term index
        bond_column_names = ["atom1", "atom2", "term_index"]

        # Create corresponding data. The first row specifies that atom0 is bonded
        # to atom 1 and has bond_type id 0
        bond_data = np.array([[0, 1, 0, ],
                              [1, 2, 0],
                              [2, 3, 0]])

        for num, name in enumerate(bond_column_names):
            bond_df[name] = bond_data[:, num]

        dl.add_bonds(bond_df)

        angle_df = pd.DataFrame()
        dihedral_df = pd.DataFrame()

        angle_column_names = ["atom1", "atom2", "atom3", "term_index"]
        dihedral_column_names = ["atom1", "atom2", "atom3", "atom4", "term_index"]

        angle_data = np.array([[0, 1, 2, 0, ],
                               [1, 2, 3, 0], ])

        dihedral_data = np.array([[0, 1, 2, 3, 0, ]])

        for num, name in enumerate(angle_column_names):
            angle_df[name] = angle_data[:, num]

        dl.add_angles(angle_df)

        for num, name in enumerate(dihedral_column_names):
            dihedral_df[name] = dihedral_data[:, num]

        dl.add_dihedrals(dihedral_df)

        if ff:
            dl.add_term_parameter(3, "harmonic", {'K': 62.100, 'theta0': 114}, uid=0,
                                  utype={'K': 'kcal * mol ** -1 * radian ** -2',
                                         'theta0': 'degree'})

            dl.add_term_parameter(2, "harmonic", {'K': 300.9, 'R0': 1.540}, uid=0,
                                  utype={'K': "kcal * mol **-1 * angstrom ** -2",
                                         'R0': "angstrom"})

        if nb:
            dl.add_nb_parameter(atom_type=1, nb_name="LJ",
                                nb_model="epsilon/sigma", nb_parameters={'sigma': 3.75, 'epsilon': 0.1947460018},
                                utype={'sigma': 'angstrom', 'epsilon': 'kcal * mol ** -1'})

            dl.add_nb_parameter(atom_type=2, nb_name="LJ",
                                nb_model="epsilon/sigma", nb_parameters={'sigma': 3.95, 'epsilon': 0.0914112887},
                                utype={'sigma': 'angstrom', 'epsilon': 'kcal * mol ** -1'})

            dl.set_mixing_rule('lorentz_berthelot')


        if scale:
            scaling_factors = {
                "coul": {
                    "scale12": 0.0,
                    "scale13": 0.0,
                    "scale14": 0.75,
                },

                "vdw": {
                    "scale12": 0.0,
                    "scale13": 0.0,
                    "scale14": 0.75,
                }
            }

            dl.set_nb_scaling_factors(scaling_factors)

        return dl

    yield build_butane
    dl.close()
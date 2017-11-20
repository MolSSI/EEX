# Plan for EEX Tests

This document outlines planned test systems for EEX. Initial testing will be in
Amber and LAMMPS.

[Energies](https://docs.google.com/spreadsheets/d/18-0ZRVw3UJ2XB_62Gu6vMm70HK1vX9669ctImiFHjfY/edit#gid=0)

[Detailed System Descriptions](https://drive.google.com/open?id=13c_BCpudbVSM34VeL7DCx1JszUu-zGTS8gCPPP39y_c)

[Information about nonbonded parameters](https://docs.google.com/document/d/1CzPa1n_d9lWwUEbFerinsNzL1R7NXQIQUykCKet4cJE/edit)

__To Do__ items are marked in __bold__.

We will first finish systems outlined in

## Tested Forcefields
1. TraPPE Forcefield
  * Harmonic Bond Parameters from the general Amber Force
Field (gaff)  
** Write about TraPPE FF here. Include references and links.  

2. General Amber Force Field (gaff)
** Write about gaff here. Include references and links.

## 1. TraPPE Forcefield - United Atom

### 1.1 Uncharged systems
1. Ethane
  - Used to test bond energy
  - __To Do__ - Test 1-2 exclusions  
2. Propane
 - Used to test angle energy
 - __To Do__ - Test 1-3 exclusions  
3. Butane
 - Used to test dihedral energy
 - __To Do__ - Make three other dihedral angles (test all quadrants)
 - __To Do__ - Test 1-4 scaling and exclusions.
4. Cyclopentane  
  - Used to test how software handles total dihedral and angle energy.
    - For AMBER, this system should result in a nonbonded energy of 0 because 1-2
  and 1-3 bonded terms are always excluded by leap. In LAMMPS, this allows for testing
  of the **special_bonds**
  - Used to test if 1-3 interactions (angle) take precedence over 1-4 nonbonded interactions.
    - In the cyclopentane molecule (TraPPE), total energy should be sum of bonded terms
    (E_NB = 0). In cyclopentane, the greatest bonded distance between atoms around
    the ring is two bonds (this is part of it, but I don't really understand yet.)
5. Cyclohexane  
 - Cyclohexane is used to test how codes handle the possibility of double counted 1-4 interactions in ring systems.
6. Long carbon chain
 - Used to test how cut-off is handled (intramolecular)
7. Two methane molecules at varying distances
   - Test cut-off schemes (LJ cut, cut+tail correction, cut + shift, smoothing
     function etc.)
     - TraPPE recommends a cut-off of 14 angstrom
     - Distances: 2 to 16 angstrom - 2 angstrom interval, & Rmin  
     \* This will have to be a periodic system (with box) since Amber doesn't
     allow cut-offs in non-periodic systems.
8. Mixture of lennard jones in the bulk
	- Test mixing rules
	- Test long tail corrections atoms different types
9. Benzene or some other molecule with impropers

### 1.2 Charged systems (electrostatics)
1. Methanol
 - Three beads with partial charges. 1-3 electrostatic interactions/exclusions.
2. Ethanol
 - Four beads with partial charges. 1-4 electrostatic interactions
3. tetrahydrofuran (five membered ring)
4. 1,3,5-trioxane (six membered ring)
5. Long chain
  - Intramolecular electrostatics
  - dimethoxyethane (CH3-O-(CH2)2-O-CH3). Varying C's in middle (1,n - dimethoxy___)
6. Charged "methanes" at varying distances
  - Test electrostatics
   * Periodic and non-periodic systems

### 2. gaff Forcefield - All Atom
All atom forcefield introduces handling of hydrogen atoms (SHAKE).
1. Nitrogen (N2)
2.

## 3. Bulk systems - Benchmarks from NIST
1. Bulk LJ fluid
2. SPCE water with long range electrostatic (ewald or pme)
3.

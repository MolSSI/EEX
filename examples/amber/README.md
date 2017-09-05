# Amber prmtop information

## File format
Amber prmtop (parameter/topology) files are separated into sections using flags with the format

%FLAG SECTION_NAME

The first section (TITLE) is one line.  From documentation - "While the title serves a primarily cosmetic purpose, this section must be present."

The second section (POINTERS) contains information about parameters in the rest of the sections.

%FLAG POINTERS
%FORMAT(10I8)
NATOM  NTYPES NBONH  MBONA  NTHETH MTHETA NPHIH  MPHIA  NHPARM NPARM
NNB    NRES   NBONA  NTHETA NPHIA  NUMBND NUMANG NPTRA  NATYP  NPHB
IFPERT NBPER  NGPER  NDPER  MBPER  MGPER  MDPER  IFBOX  NMXRS  IFCAP
NUMEXTRA NCOPY

More information can be found at the resources below:

* Reading and modifying Amber parameter files - http://alma.karlov.mff.cuni.cz/bio/99_Studenti/00_Dalsi/ParamFit/2013_ParamFit_AmberTools13.pdf

* Information on prmtop format from ambermd.org - http://ambermd.org/prmtop.pdf


## FF Parameter Units

This section lists the type of unit associated with flags which specify the force
field parameters (or the "energy expression"). Other flags not listed specify topology

Flag                     | Symbol            | Units                   | Special notes
 -----------------       | :--------: | :---------------------: | :-----------
 BOND_FORCE_CONSTANT     | k          | kcal/mol/Angstrom$^2    | (1/2)k(r-r_eq)
 BOND_EQUIL_VALUE        | r_eq       | Angstrom                | (1/2)k(r-r_eq)
 ANGLE_FORCE_CONSTANT    | k_theta    | kcal/mol/radian$^2$     | (1/2)k_theta(theta-theta_eq)
 ANGLE_EQUIL_VALUE       | theta_eq   | radians                 |
 DIHEDRAL_FORCE_CONSTANT | k_tor      | kcal/mol                |
 DIHEDRAL_PERIODICITY    | n          |                         |
 DIHEDRAL_PHASE          |            |                         |
 LENNARD_JONES_ACOEF     | A_ij       | Angstrom^12 kcal/mol    | Used in prmtop file only
 LENNARD_JONES_BCOEF     | B_ij       | Angstrom^6 kcal/mol     | Used in prmtop file only
 SCEE_SCALE_FACTOR       |            |                         | Scaling factor for 1-4 electrostatics
SCNB_SCALE_FACTOR       |            |                         | Scaling factor for 1-4 onbonded interactions

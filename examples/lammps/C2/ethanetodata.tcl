#!/usr/bin/tclsh

# topotools example:
# generate topology data from coordinate data
# build a box of 8x8x8 = 512 molecules.
##############################################

# explicitly load topotools and pbctools packages since
# they are not automatically requred in text mode and
# abort if their version numbers are insufficient.
if {[catch {package require topotools 1.1} ver]} {
   vmdcon -error "$ver. This script requires at least TopoTools v1.1. Exiting..."
   quit
}

if {[catch {package require pbctools 2.3} ver]} {
   vmdcon -error "$ver. This script requires at least pbctools v2.3. Exiting..."
   quit
}

set fname ethane.xyz
# check for presence of coordinate file
if {! [file exists $fname]} {
   vmdcon -error "Required file '$fname' not available. Exiting..."
   quit
}
  
# load coordinates and use automatically computed bonds
mol new $fname autobonds yes waitfor all

# create separate selections for hydrogens and carbons
# and one for all atoms.
set sel [atomselect top all]
$sel set type CH3
$sel set mass 15.03452
$sel set charge 0.000

# with a proper .pdb file, VMD will have already
# determined the proper element definitions, so
# recomputing the bonds will be hardly necessary.
# we still need to assign bond types, though.
topo retypebonds 
vmdcon -info "assigned [topo numbondtypes] bond types to [topo numbonds] bonds:"
vmdcon -info "bondtypes: [topo bondtypenames]"

# now derive angle and dihedral definitions from bond topology.
# every two bonds that share an atom yield an angle.
# every two bonds that share a bond yield a dihedral.
topo guessangles
vmdcon -info "assigned [topo numangletypes] angle types to [topo numangles] angles:"
vmdcon -info "angletypes: [topo angletypenames]"
topo guessdihedrals
vmdcon -info "assigned [topo numdihedraltypes] dihedral types to [topo numdihedrals] dihedrals:"
vmdcon -info "dihedraltypes: [topo dihedraltypenames]"

# now let VMD reanalyze the molecular structure
# this is needed to detect fragments/molecules
# after we have recomputed the bonds
mol reanalyze top

## now set box dimensions from the min/max corners in order 
## to fit the molecule  considering its vdw atom radii.
#set minmax [measure minmax $sel -withradii]
## we need to increase the box by 10% to get a reasonable density.
#set box [vecscale 1.1 [vecsub [lindex $minmax 1] [lindex $minmax 0]]]
#pbc set $box
#vmdcon -info "box size: $box"
## and recenter the coordinates around the center of mass
#set center [measure center $sel weight none]
#$sel moveby [vecscale -1.0 $center]
#vmdcon -info "moved center from $center to [measure center $sel weight none]"

# now set box dimensions and write out the result as 
# a lammps data file.
pbc set {10.0 10.0 10.0 90.0 90.0 90.0}

# we use a high-level tool from to multiply the system.
TopoTools::replicatemol top 1 1 1

# and write out the result as a lammps data file.
topo writelammpsdata data.ethane full

## for easier testing and visualization, we
## also write out copies in .pdb and .psf format.
#animate write pdb 64xmethane.pdb
#animate write psf 64xmethane.psf
#
## done. now exit vmd
quit

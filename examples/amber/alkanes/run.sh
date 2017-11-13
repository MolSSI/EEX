#!/bin/bash 

# Remove old energy files
rm energies.txt

# System names

declare -a arr=("ethane" "propane" "butane" "butane_conformation2")

# Start energy calculation file
echo "Single point energy calculations" > energies.txt

for i in "${arr[@]}"
do

	cat > leap.in <<EOL
	loadamberparams frcmod_trappe.$(echo ${i})
	sys = loadmol2 $(echo ${i}).mol2
	saveamberparm sys trappe_$(echo ${i})_single_molecule.prmtop trappe_$(echo ${i})_single_molecule.inpcrd
	quit
	
EOL
	
	tleap -f leap.in
	
	cpptraj trappe_$(echo ${i})_single_molecule.prmtop<<!
	trajin trappe_$(echo ${i})_single_molecule.inpcrd
	trajout trappe_$(echo ${i})_single_molecule.nc
	go
!

	
	sander -O -i single_point.in -c trappe_$(echo ${i})_single_molecule.inpcrd -p trappe_$(echo ${i})_single_molecule.prmtop -ref trappe_$(echo ${i})_single_molecule.inpcrd -o single_point_$(echo ${i}).out -r out.rst7 -y trappe_$(echo ${i})_single_molecule.nc
	
    echo "$i" >> energies.txt
    awk '/NSTEP/, /minimization/' single_point_$i.out >> energies.txt
    echo "" >> energies.txt
done

rm *.nc
rm leap.in
rm *.out


#!/bin/bash

# Remove old energy files
rm energies*.txt
rm energies.csv

# System names

declare -a arr=("methanol" "ethanol" )

# Start energy calculation file
# echo "Single point energy calculations" > energies.txt
# echo "${arr[@]}" >> energies.txt

echo "molecule,bond,angle,dihedral,vdwaals,electrostatic,hbond,1-4_VDW,1-4_EEL,restraint," > energies.csv

for i in "${arr[@]}"
do

	cpptraj trappe_$(echo ${i})_single_molecule.prmtop<<!
	trajin trappe_$(echo ${i})_single_molecule.inpcrd
	trajout trappe_$(echo ${i})_single_molecule.nc
	go
!


	sander -O -i single_point.in -c trappe_$(echo ${i})_single_molecule.inpcrd -p trappe_$(echo ${i})_single_molecule.prmtop -ref trappe_$(echo ${i})_single_molecule.inpcrd -o single_point_$(echo ${i}).out -r out.rst7 -y trappe_$(echo ${i})_single_molecule.nc

    #echo "$i" >> energies.txt
    awk '/ BOND /, / RESTRAINT /' single_point_$i.out >> energies_$i.txt

		# Replace spaces in 1-4 VDW and 1-4 EEL with underscores for easier text
		# processing

		perl -i -pe's/1-4 VDW/1-4_VDW/g' energies_$i.txt
		perl -i -pe's/1-4 EEL/1-4_EEL/g' energies_$i.txt
    echo "" >> energies_$i.txt

		# Process energies.txt to CSV
		echo -n "trappe_$(echo ${i})_single_molecule," >> energies.csv
		awk -v ORS="," 'FNR == 1 {print $3}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 1 {print $6}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 1 {print $9}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 2 {print $3}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 2 {print $6}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 2 {print $9}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 3 {print $3}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 3 {print $6}' energies_$i.txt >> energies.csv
		awk -v ORS="," 'FNR == 3 {print $9}' energies_$i.txt >> energies.csv
		echo "" >> energies.csv
done

rm *.nc
rm *.out
rm *.txt

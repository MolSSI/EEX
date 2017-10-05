cpptraj trappe_ethane_single_molecule.prmtop<<!
trajin trappe_ethane_single_molecule.inpcrd
trajout trappe_ethane_single_molecule.nc
go
!

cpptraj trappe_propane_single_molecule.prmtop<<!
trajin trappe_propane_single_molecule.inpcrd
trajout trappe_propane_single_molecule.nc
go
!

cpptraj trappe_butane_single_molecule.prmtop<<!
trajin trappe_butane_single_molecule.inpcrd
trajout trappe_butane_single_molecule.nc
go
!

sander -O -i single_point.in -c trappe_ethane_single_molecule.inpcrd -p trappe_ethane_single_molecule.prmtop -ref trappe_ethane_single_molecule.inpcrd -o single_point_ethane.out -r out.rst7 -y trappe_ethane_single_molecule.nc

sander -O -i single_point.in -c trappe_propane_single_molecule.inpcrd -p trappe_propane_single_molecule.prmtop -ref trappe_propane_single_molecule.inpcrd -o single_point_propane.out -r out.rst7 -y trappe_propane_single_molecule.nc

sander -O -i single_point.in -c trappe_butane_single_molecule.inpcrd -p trappe_butane_single_molecule.prmtop -ref trappe_butane_single_molecule.inpcrd -o single_point_butane.out -r out.rst7 -y trappe_butane_single_molecule.nc

#sander -O -i min.in -c trappe_ethane_single_molecule.inpcrd -p trappe_ethane_single_molecule.prmtop -ref trappe_ethane_single_molecule.inpcrd -o single_point_ethane.out -r trappe_ethane_single_molecule.rst7

#sander -O -i min.in -c trappe_propane_single_molecule.inpcrd -p trappe_propane_single_molecule.prmtop -ref trappe_propane_single_molecule.inpcrd -o single_point_propane.out -r trappe_propane_single_molecule.rst7

#sander -O -i min.in -c trappe_butane_single_molecule.inpcrd -p trappe_butane_single_molecule.prmtop -ref trappe_butane_single_molecule.inpcrd -o single_point_butane.out -r trappe_butane_single_molecule.rst7

# Extract energies and clean
declare -a arr=("ethane" "propane" "butane")

echo "Single point energy calculations" > energies.txt

for i in "${arr[@]}"
do
    echo "$i" >> energies.txt
    awk '/NSTEP/, /minimization/' single_point_$i.out >> energies.txt
    echo "" >> energies.txt
done

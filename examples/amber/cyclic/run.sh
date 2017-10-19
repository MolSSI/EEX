cpptraj trappe_cyclopentane_single_molecule.prmtop<<!
trajin trappe_cyclopentane_single_molecule.inpcrd
trajout trappe_cyclopentane_single_molecule.nc
go
!

cpptraj trappe_cyclohexane_single_molecule.prmtop<<!
trajin trappe_cyclohexane_single_molecule.inpcrd
trajout trappe_cyclohexane_single_molecule.nc
go
!

sander -O -i single_point.in -c trappe_cyclopentane_single_molecule.inpcrd -p trappe_cyclopentane_single_molecule.prmtop -ref trappe_cyclopentane_single_molecule.inpcrd -o single_point_cyclopentane.out -r out.rst7 -y trappe_cyclopentane_single_molecule.nc

sander -O -i single_point.in -c trappe_cyclohexane_single_molecule.inpcrd -p trappe_cyclohexane_single_molecule.prmtop -ref trappe_cyclohexane_single_molecule.inpcrd -o single_point_cyclohexane.out -r out.rst7 -y trappe_cyclohexane_single_molecule.nc


# Extract energies and clean
declare -a arr=("cyclopentane" "cyclohexane" )

echo "Single point energy calculations" > energies.txt

for i in "${arr[@]}"
do
    echo "$i" >> energies.txt
    awk '/NSTEP/, /minimization/' single_point_$i.out >> energies.txt
    echo "" >> energies.txt
done

rm *.nc

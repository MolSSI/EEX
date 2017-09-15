# Units based LAMMPS unit styles (http://lammps.sandia.gov/doc/units.html)


_units_style = {
    "lj": {},
    "real": {
        "mass": "gram  mol^-1",
        "distance" : "Angstrom",
        "time" : "femtosecond",
        "energy" : "kcal / mol",
        "velocity" : "Angstrom / femtosecond",
        "force" : "kcal  mol^-1 Angstrom^-1",
        "torque" : "kcal  mol^-1",
        "temperature" : "Kelvin",
        "pressure" : "atmosphere",
        "dynamic viscosity" : "poise",
        "charge" : "e",
        "dipole" : "e  Angstrom",
        "electric field" : "volt  Angstrom^-1",
        "density" : "gram  cm^-dim"

    },

    "metal": {

        "mass" : "gram mol^-1",
        "distance" : "Angstrom",
        "time" : "picosecond",
        "energy" : "eV",
        "velocity" : "Angstrom  picosecond^-1",
        "force" : "eV Angstrom^-1",
        "torque" : "eV",
        "temperature" : "Kelvin",
        "pressure" : "bar",
        "dynamic viscosity" : "poise",
        "charge" : "e",
        "dipole" : "e * Angstrom",
        "electric field" : "volt  Angstrom^-1",
        "density" : "gram  cm^-dim",
    },

    "si": {

        "mass" : "kilogram",
        "distance" : "meter",
        "time" : "second",
        "energy" : "Joule",
        "velocity" : "meter second^-1",
        "force" : "Newtons",
        "torque" : "Newton  meter",
        "temperature" : "Kelvin",
        "pressure" : "Pascal",
        "dynamic viscosity" : "Pascal  second",
        "charge" : "Coulomb",
        "dipole" : "Coulomb meter",
        "electric field" : "volt meter^-1",
        "density" : "kilogram meter^-dim",
    },


    "cgs": {
        "mass" : "gram",
        "distance" : "centimeters",
        "time" : "second",
        "energy" : "ergs",
        "velocity" : "centimeters second^-1",
        "force" : "dyne",
        "torque" : "dyne - centimeters",
        "temperature" : "Kelvin",
        "pressure" : "dyne cm^-2",
        "dynamic viscosity" : "poise",
        "charge" : "statcoulomb",
        "dipole" : "statcoulomb centimeter",
        "electric field" : "statvolt cm^-1",
        "density" : "gram cm^-dim,
    },

    "electron": {

        "mass" : "AMU",
        "distance" : "Bohr",
        "time" : "femtosecond",
        "energy" : "Hartree",
        "velocity" : "Bohr (atomic time unit)^-1",
        "force" : "Hartree Bohr^-1",
        "temperature" : "Kelvin",
        "pressure" : "Pascal",
        "charge" : "e",
        "dipole moment" : "Debye",
        "electric field" : "volt cm^-1",
},

    "micro": {

        "mass" : "picogram",
        "distance" : "micrometer",
        "time" : "microsecond",
        "energy" : "picogram  micrometer^2 microsecond^-2",
        "velocity" : "micrometers microsecond^-1",
        "force" : "picogram micrometer microsecond^-2",
        "torque" : "picogram micrometer^2 microsecond^-2",
        "temperature" : "Kelvin",
        "pressure" : "picogram  micrometer^-1 microsecond^-2",
        "dynamic viscosity" : "picogram micrometer^-1 microsecond^-1",
        "charge" : "picocoulomb",
        "dipole" : "picocoulomb micrometer",
        "electric field" : "volt micrometer^-1",
        "density" : "picogram micrometer^-dim"

    },

    "nano": {
        "mass" : "attogram",
        "distance" : "nanometer",
        "time" : "nanosecond",
        "energy" : "attogram nanometer^2 nanosecond^-2",
        "velocity" : "nanometers nanosecond^-1",
        "force" : "attogram nanometer nanosecond^-2",
        "torque" : "attogram nanometer^2 nanosecond^-2",
        "temperature" : "Kelvin",
        "pressure" : "attogram nanometer^-1 nanosecond^-2",
        "dynamic viscosity" : "attogram nanometer^-1 nanosecond^-1",
        "charge" : "e",
        "dipole" : "e nanometer",
        "electric field" : "volt nanometer^-1",

    },

}
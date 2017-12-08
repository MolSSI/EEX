import eex
import eex_find_files

def build_dl(program, test_directory, molecule):
    if program.lower() == "amber":
        file_name = "%s.prmtop" % molecule
        fname = eex_find_files.get_example_filename("amber", test_directory, file_name)
        dl = eex.datalayer.DataLayer("test_amber")
        data = eex.translators.amber.read_amber_file(dl, fname)
        return data, dl
    elif program.lower() == "lammps":
        file_name = molecule
        fname = eex_find_files.get_example_filename("lammps", test_directory, file_name)
        dl = eex.datalayer.DataLayer("test_lammps")
        data = eex.translators.lammps.read_lammps_file(dl, fname)
        return data, dl
    else:
        raise KeyError("Program %s not understood" % program)

import numpy as np
import pandas as pd
import os
from glob import glob
import peakutils
import yaml


def main():
    """ Main driver function for the script.

        The script will find a input.yml to use for
        all of the parameters, and then look for
        .xyz files to spawn the calculations.
    """
    print("Reading in settings")

    with open("input.yml") as read_file:
        calc_params = yaml.load(read_file)

    with open("template.com") as read_file:
        template = read_file.read()

    print("Read in the following:")
    print("----------------------")
    for key, value in calc_params.items():
        print(key, value)

    print("----------------------")

    xyzfiles = glob("*.xyz")
    print("Found the following .xyz")
    print(xyzfiles)

    if os.path.isdir("calcs") is False:
        os.mkdir("calcs")

    for file in xyzfiles:
        natoms, comment, coords = read_xyz(file)
        # Automatically make doublets based on the number
        # of hydrogens - usually even number of hydrogens is
        # closed-shell.
        hcount = len([line[0] for line in coords if line[0] == "H"])
        if hcount % 2 != 0:
            calc_params["multiplicities"] = [2]
        elif hcount % 2 == 0:
            calc_params["multiplicities"] = [1]
        input_builder(
            template,
            coords,
            calc_params,
            comment)


def read_xyz(filepath):
    """ Parses an .xyz file, and returns the
        number of atoms and the comment line.
    """
    with open(filepath) as read_file:
        lines = read_file.readlines()
        natoms = int(lines[0])
        # Take the comment line, and get rid of newline characters
        comment = lines[1].replace("\n", "")
        comment = comment.replace(" ", "")
        coordinates = lines[2:]
        return natoms, comment, coordinates


def input_builder(template, coords, params, comment):
    """ Function for building a Gaussian input file.
        Takes a template input file, 
    """
    coords = "".join(coords)
    for multiplicity in params["multiplicities"]:
        # Flatten the coordinates into space delimited
        if params["ts"]:
            opt = "Opt=(VeryTight,CalcAll,TS,MaxCycles=100)"
        else:
            opt = "Opt=(MaxCycles=100,VeryTight)"
        # Package parameters together
        filename = comment + "-" + str(multiplicity) + ".com"
        name = filename.replace(".com", "")
        params.update({
            "multi": str(multiplicity),
            "comment": comment,
            "coords": coords,
            "opt": opt,
            "name": name
            })
        # Write the input file to disk
        with open("calcs/" + filename, "w+") as write_file:
            write_file.write(template.format_map(params))
        with open("g16.sh") as read_file:
            pbs_template = read_file.read()
        with open("calcs/g16.sh", "w+") as write_file:
            write_file.write(
                pbs_template.format_map({"name": name})
                )
        os.chdir("calcs")
        os.system("qsub g16.sh")
        os.chdir("..")


if __name__ == "__main__":
    main()

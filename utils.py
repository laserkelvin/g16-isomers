
"""
    Contains all of the auxillary functions for CFOURviewer; i.e. file I/O,
    storage to HDF5, copying and pasting, etc.

    The scope of HDF5 will be to store the parsed data, as well as the full
    output file as a string.

    Settings will be stored in a dot folder in the user's home directory;
    this includes templates for the PBS script as well as CFOURviewer settings.
"""

import os
import shutil
import datetime
import re
from dataclasses import dataclass

import numpy as np
from subprocess import Popen, run, PIPE
from glob import glob
from itertools import product
import h5py
import yaml
import joblib
import periodictable


def generate_folder():
    # Function to generate the next folder in the ID chain.
    # Returns the next ID number in the chain as an integer.
    settings = read_settings()
    dir_list = glob(settings["calc_dir"] + "/*")
    filtered = list()
    for dir in dir_list:
        # This takes only folder names that are numeric
        try:
            filtered.append(int(dir))
        except TypeError:
            pass
    next_ID = max(filtered) + 1
    os.mkdir(settings["calc_dir"] + str(next_ID))
    return next_ID


def read_settings():
    # Wrapper for the read_yaml function, specifically to call up the
    # settings file.
    location = os.path.expanduser("~") + "/.cfourviewer/settings.yml"
    return read_yaml(location)

"""
    File I/O

    Includes YAML and HDF5 functions

    HDF5 system is organized into IDs - an ID can contain one or several
    calculations, and the attributes of an ID group are metadata regarding
    the calculation batch, i.e. a SMILES code to identify the molecule.

    Each calculation is then stored as datasets within this group, and the
    parsed results of the calculation.
"""

def write_yaml(yaml_path, contents):
    # Function for writing dictionary to YAML file
    with open(yaml_path, "w+") as write_file:
        yaml.dump(contents, write_file, default_flow_style=False)


def read_yaml(yaml_path):
    # Function for reading in a YAML file
    with open(yaml_path) as read_file:
        return yaml.load(read_file)


def xyz2smi(filepath):
    # Calls external obabel executable to convert an input xyz file
    # Returns a SMILES string for identification
    proc = Popen(
        ["obabel", "-ixyz", filepath, "-osmi"],
        stdout=PIPE
    )
    output = proc.communicate()[0].decode()
    smi = output.split("\t")[0]
    return smi


def log2xyz(filepath):
    # Calls external obabel executable to convert g09 logfile
    # Dumps the xyz to disk
    command = ["obabel", "-ig03", filepath, "-oxyz"]
    with Popen(command, stdout=PIPE) as babel_proc:
        output = babel_proc.communicate()[0].decode()
    with open("xyz", "w+") as write_file:
        write_file.write(output)
    return output


def read_xyz(filepath):
    """read_xyz

    :param filepath:
    """
    natoms_re = re.compile(r"^\d*\n", re.M)
    coords_re = re.compile(r"^[A-Z]\s*[-]?\d.\d*\s*[-]?\d.\d*\s*[-]?\d.\d*", re.M)
    with open(filepath) as read_file:
        contents = read_file.read()
        natoms = int(natoms_re.findall(contents)[0])
        coords = coords_re.findall(contents)
        coords = [coord.split() for coord in coords]
    mol_dict = {
            "natoms": natoms,
            "comment": "",
            "coords": coords
            }
    return mol_dict


def xyz2str(xyz_list):
    """ Convert xyz from a list format to a
        formatted string
    """
    # flatten each row into string
    xyz_str = [" ".join(row) for row in xyz_list]
    # flatten rows into a single string
    xyz_str = "\n".join(xyz_str)
    return xyz_str


"""
    Miscellaneous functions
"""

def combine_method_basis(methods, bases):
    """ Function to return an iterator over all possible
        combinations of methods and bases.
    """
    return product(methods, bases)


def run_g16():
    # Wrap g09 executable
    command = ["g16", "calc.com"]
    with open("calc.log", "w+") as logfile:
        with open("calc.com", "r") as comfile:
            with Popen(command, stdin=comfile, stdout=logfile) as g16_proc:
                g16_proc.wait()


def check_calc(method, basis):
    """ Function to see if calculation has already been
        performed. Returns a boolean indicating False
        for not done, and True for completed calcs.
    """
    # If there is no logfile then it can't have been done
    if os.path.isfile("done") is False:
        return False
    else:
        # Open the log file to check
        with open("done") as read_file:
            # If the method/basis is logged, then the
            # calculation has been completed
            text = read_file.read()
            if "/".join([method, basis]) in text:
                return True
            else:
                return False
    

def smi2xyz(smi_file):
    # Generate structures and optimize with UFF into a big SDF file
    bab_cmd = f"obgen {smi_file} -ff UFF -n 200"
    with open("structures/full.sdf", "w+") as write_file:
        babel_proc = run(bab_cmd, shell=True, stdout=write_file)
        babel_proc.wait()
    os.chdir("structures")
    convert_cmd = "obabel -isdf full.sdf -O geom.xyz -m"
    babel_proc = Popen(convert_cmd, shell=True, stdout=PIPE)
    babel_proc.wait()
    os.chdir("..")


def save_obj(obj, filepath, **kwargs):
    """
    Save an object to disk.
    """
    settings = {
        "compress": ("gzip", 9)
    }
    settings.update(kwargs)
    joblib.dump(obj, filepath, **kwargs)


@dataclass
class Molecule:
    A: float = 0.0
    B: float = 0.0
    C: float = 0.0
    success: bool = False
    u_A: float = 0.0
    u_B: float = 0.0
    u_C: float = 0.0
    formula: str = ""
    smi: str = ""
    point_group: str = "C1"
    method: str = ""
    basis: str = ""
    charge: int = 0
    multi: int = 1
    kappa: float = 0.0
    DJ: float = 0.0
    DJK: float = 0.0
    DK: float = 0.0
    delJ: float = 0.0
    delK: float = 0.0
    Iaa: float = 0.0
    Ibb: float = 0.0
    Icc: float = 0.0
    defect: float = 0.0
    coords: str = ""
    zpe: float = 0.0
    Etot: float = 0.0
    harm_freq: str = ""
    harm_int: str = ""
    opt_delta: float = 0.0
    filename: str = ""

    def __eq__(self, other, thres=1e-3):
        check = all(
            [
                np.abs(self.A - other.A) <= thres,
                np.abs(self.B - other.B) <= thres,
                np.abs(self.C - other.C) <= thres
                ]
            )
        return check


def parse_g16(filepath):
    """parse_g16

    :param filepath: Path to the logfile
    """
    data = dict()
    harm_freq = list()
    harm_int = list()
    filename = filepath.split("/")[-1].split(".")[0]
    with open(filepath) as read_file:
        lines = read_file.readlines()
        for index, line in enumerate(lines):
            if "Rotational constants (MHZ)" in line:
                rot_con = lines[index + 1].split()
                rot_con = [float(value) for value in rot_con]
                A, B, C = rot_con
                data["A"] = A
                data["B"] = B
                data["C"] = C
            if "Dipole moment (Debye)" in line:
                dipoles = lines[index + 1].split()[:3]
                dipoles = [float(value) for value in dipoles]
                u_A, u_B, u_C = dipoles
                data["u_A"] = u_A
                data["u_B"] = u_B
                data["u_C"] = u_C
            if "Full point group" in line:
                data["point_group"] = line.split()[3]
            if "Stationary point found" in line:
                data["success"] = True
            if line.startswith(" # "):
                calc = line.split()[1].split("/")
                try:
                    method, basis = calc
                    data["basis"] = basis
                except ValueError:
                    # This is for composite schemes
                    method = calc
                data["method"] = method
            if "Multiplicity" in line:
                split_line = line.split()
                data["charge"] = int(split_line[2])
                data["multi"] = int(split_line[-1])
            if "Vibro-Rot alpha Matrix" in line:
                alpha_flag = True
                alpha_lines = lines[index + 3:]
                alpha_mat = list()
                alpha_index = 0
                while alpha_flag is True:
                    current_line = alpha_lines[alpha_index]
                    if current_line.startswith("Q("):
                        alpha = alpha_lines[alpha_index].split()[2:]
                        alpha = [float(value) for value in alpha]
                        alpha_mat.append(alpha)
                        alpha_index += 1
                    else:
                        alpha_flag = False
            if "Asymm. param." in line:
                data["kappa"] = float(line.split()[-1])
            if "DELTA J  :" in line:
                data["DJ"] = float(line.replace("D", "E").split()[-1])
            if "DELTA JK :" in line:
                data["DJK"] = float(line.replace("D", "E").split()[-1])
            if "DELTA K  :" in line:
                data["DK"] = float(line.replace("D", "E").split()[-1])
            if "delta J  :" in line:
                data["delJ"] = float(line.replace("D", "E").split()[-1])
            if "delta K  :" in line:
                data["delK"] = float(line.replace("D", "E").split()[-1])
            if "Iaa" in line:
                split_line = line.replace("D", "E").split()
                data["Iaa"] = float(split_line[2])
                data["Ibb"] = float(split_line[4])
                data["Icc"] = float(split_line[-1])
                data["defect"] = data["Icc"] - data["Iaa"] - data["Ibb"]
            if "Principal axis orientation" in line:
                coord_lines = lines[index + 5:]
                coord_flag = True
                coord_mat = list()
                coord_index = 0
                while coord_flag is True:
                    current_line = coord_lines[coord_index]
                    if "------" in current_line:
                        coord_flag = False
                    else:
                        coords = current_line.split()[1:]
                        coords = [float(value) for value in coords]
                        coord_mat.append(coords)
                        coord_index += 1
                data["coords"] = np.array(coord_mat)
            if "Zero-point correction" in line:
                data["zpe"] = float(line.split()[2])
            if "Sum of electronic and zero-point" in line:
                data["Etot"] = float(line.split()[-1])
            if "Frequencies --" in line:
                freq = line.split()[2:]
                freq = [float(value) for value in freq]
                harm_freq.extend(freq)
            if "IR Inten" in line:
                inten = line.split()[3:]
                inten = [float(value) for value in inten]
                harm_int.extend(inten)
            if "Predicted change in Energy=" in line:
                data["opt_delta"] = float(line.replace("D", "E").split("=")[-1])
    if "coords" in data:
        atom_dict = dict()
        for coord in data["coords"]:
            element = periodictable.elements[coord[0]]
            if element not in atom_dict:
                atom_dict[element] = 1
            else:
                atom_dict[element] += 1
        molecule_string = "".join(["{}{}".format(key, value) for key, value in atom_dict.items()])
        data["formula"] = molecule_string
    data["filename"] = filename
    data["harm_freq"] = harm_freq
    data["harm_int"] = harm_int
    result_obj = Molecule(**data)
    return result_obj


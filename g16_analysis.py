
from glob import glob

import numpy as np
import pandas as pd


def parse_log(filepath):
    """ Function for parsing a G16 logfile.
    """
    frequencies = list()
    filename = filepath.split("/")[-1].split(".log")[0].split("-")
    isomer = filename[0]
    multi = filename[1]
    ZPE = 0.
    total = 0.
    G3 = 0.
    success = False
    with open(filepath) as read_file:
        for line in read_file:
            if "Zero-point correction=" in line:
                ZPE = float(line.split()[2])
            if "Sum of electronic and zero-point" in line:
                total = float(line.split()[-1])
            if "Frequencies --" in line:
                current = line.split()[2:]
                frequencies.extend([float(freq) for freq in current])
            if "Normal termination" in line:
                success = True
            if "G3(0 K)" in line:
                G3 = float(line.split()[2])
    data = {
        "Isomer": isomer,
        "Multiplicity": multi,
        "success": success,
        "ZPE": ZPE,
        "G3-0K": G3,
        "E+ZPE": total,
        "Frequencies": frequencies
        }
    return data


def main(directory="calcs"):
    conversion = 219474.6 / 83.59
    logfiles = glob(os.path.join(directory, "*.log"))
    data = [parse_log(file) for file in logfiles]
    df = pd.DataFrame(data)
    df["E"] = df["E+ZPE"] - df["ZPE"]
    df["Relative-0K"] = conversion * (df["E+ZPE"] - df["E+ZPE"].min())
    df["Relative-E"] = conversion * (df["E"] - df["E"].min())
    df["Relative-G3"] = conversion * (df["G3-0K"] - df["G3-0K"].min())
    df.to_csv("collected_results.csv")


if __name__ == "__main__":
    main()

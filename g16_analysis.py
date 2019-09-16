#!/usr/bin/env python

from glob import glob
from pathlib import Path
import tarfile
import os

import pandas as pd
from tqdm.autonotebook import tqdm

from utils import parse_g16, save_obj


def main():
    data = list()
    ignore = ["harm_freq", "harm_int"]
    molecules = list()
    calc_path = Path("calcs")
    if calc_path.exists() is False:
        raise Exception("No calculation folder!")
    for file in tqdm(calc_path.rglob("*.log")):
        molecule = parse_g16(file)
        if molecule not in molecules:
            molecules.append(molecule)
    for molecule in molecules:
        data.append(
            {key: value for key, value in molecule.__dict__.items() if key not in ignore}
            )
    mol_df = pd.DataFrame(data)
    mol_df.sort_values(["Etot"], ascending=True, inplace=True)
    indexes = mol_df.index
    # Resort the molecule list with ascending order energy
    molecules = [molecules[index] for index in indexes]
    mol_df.loc[:, "relative"] = (mol_df["Etot"] - mol_df["Etot"].min()) * (219474.3 / 83.59)
    mol_df.to_pickle("mol_dataframe.pkl")
    success_df = mol_df.loc[mol_df["success"] == True]
    success_df.to_pickle("mol_dataframe-success.pkl")
    save_obj(molecules, "molecules.pkl")

    with open("result_template.txt", "r") as read_file:
        template = read_file.read()

    with open("result_summary.txt", "w+") as write_file:
        for index, molecule in enumerate(molecules):
            molecule.rel_energy = mol_df.iloc[index]["relative"]
            write_file.write(template.format(**molecule.__dict__))
    
    with tarfile.open("compressed_calcs.tar.gz", "w:gz") as tar_file:
        for file in calc_path.rglob("*.log"):
            tar_file.add(file)
        tar_file.close()

if __name__ == "__main__":
    main()


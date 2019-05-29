
from glob import glob
import os

import pandas as pd
from tqdm.autonotebook import tqdm

from utils import parse_g16, save_obj


def main():
    data = list()
    ignore = ["coords", "harm_freq", "harm_int"]
    molecules = [parse_g16(file) for file in tqdm(glob("calcs/*.log"))]
    for molecule in molecules:
        data.append(
            {key: value for key, value in molecule.__dict__.items() if key not in ignore}
            )
    mol_df = pd.DataFrame(data)
    mol_df.to_pickle("mol_dataframe.pkl")
    success_df = mol_df.loc[mol_df["success"] == True]
    success_df.to_pickle("mol_dataframe-success.pkl")
    save_obj(molecules, "molecules.pkl")

    with open("result_template.txt", "r") as read_file:
        template = read_file.read()

    for molecule in molecules:
        print(template.format(**molecule.__dict__))

if __name__ == "__main__":
    main()


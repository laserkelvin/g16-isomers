
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
            {key: value for key, value in molecule.__dict__ if key in ignore}
            )
    mol_df = pd.DataFrame(data)
    mol_df.to_pickle("mol_dataframe.pkl")
    save_obj(molecules, "molecules.pkl")


if __name__ == "__main__":
    main()


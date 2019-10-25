#!/usr/bin/env python

from itertools import product
from concurrent.futures import ProcessPoolExecutor

import click
import pandas as pd

import gen_pubchem
import g16_isomers
import utils


def generate_stochiometry(atom, n_min=2, n_max=8, step=1):
    """
    Function to generate combinations of atom numbers for a given
    range of possible atom numbers.
    """
    if n_min != n_max:
        combinations = [
            f"{atom}{index}" if index != 1 else f"{atom}"
            for index in range(n_min, n_max + 1, step)
        ]
    else:
        # In the instance we're locking the number of atoms
        if n_min == 1:
            combinations = [f"{atom}"]
        else:
            combinations = [f"{atom}{n_min}"]
    return combinations


def generate_formulas(chem_dict, filepath="formulas.txt"):
    atom_numbers = [
        generate_stochiometry(key, **value) for key, value in chem_dict.items()
    ]
    combinations = list()
    with open(filepath, "w+") as write_file:
        for combination in product(*atom_numbers):
            write_file.write("".join(combination) + "\n")
            combinations.append(combination)
    return combinations


@click.command()
@click.argument("yml_path")
@click.option(
    "--filepath", default="formulas.txt", help="File to save the formulas to."
)
@click.option(
    "-s", "--skip-formula", "skip", default=False, help="Skip formulae generation."
)
@click.option(
    "-n", "--nprocess", "nprocess", default=4, help="Number of concurrent processes."
)
def main(yml_path, filepath, skip, nprocess):
    if not skip:
        chem_dict = utils.read_yaml(yml_path)
        print("Generating combinations.")
        combinations = generate_formulas(chem_dict, filepath)
    else:
        with open(filepath) as read_file:
            combinations = read_file.readlines()
        combinations = [combo.replace("\n", "") for combo in combinations]
    # Run a pool of processes to speed up the querying
    with ProcessPoolExecutor(max_workers=nprocess) as pool:
        full_batch = list(
            pool.map(
                gen_pubchem.exhaust_query, ["".join(combo) for combo in combinations]
            )
        )
    # Combine all the results
    full_df = pd.concat(full_batch)
    full_df.to_csv("full_pull.csv")
    print("Generating XYZ files.")
    gen_pubchem.df2xyz(full_df)


if __name__ == "__main__":
    main()

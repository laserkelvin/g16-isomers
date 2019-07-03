#!/usr/bin/env python


import os

import pubchempy
import pandas as pd
import click

import utils


def query_pubchem(formula, max_count=100):
    df = pubchempy.get_compounds(
            formula,
            "formula",
            listkey_count=max_count,
            as_dataframe=True
            )
    return df


def exhaust_query(formula, max_count=100):
    """
    This function will persist in querying until no more
    molecules are known in the database.
    """
    batch = list()
    pubchem_df = query_pubchem(formula, max_count)
    batch.append(pubchem_df)
    start = len(pubchem_df)
    while len(pubchem_df) == max_count:
        if start >= 500:
            break
        pubchem_df = query_pubchem(formula, max_count)
        start += len(pubchem_df)
        batch.append(pubchem_df)
    # Combine the dataframes together
    full_df = pd.concat(batch)
    return full_df
    

def df2xyz(pubchem_df):
    # Write smiles to file
    smi = pubchem_df["canonical_smiles"].values
    with open("molecules.smi", "w+") as write_file:
        for value in smi:
            write_file.write(value + "\n")
    if os.path.isdir("structures") is False:
        os.mkdir("structures")
    utils.smi2xyz("molecules.smi")



@click.command()
@click.argument("formula")
@click.option("--max_count", default=100, help="Maximum number of queries to return.")
def main(formula, max_count=100):
    df = query_pubchem(formula, max_count)
    df.to_csv("pubchem_query.csv")
    df2xyz(df)


if __name__ == "__main__":
    main()


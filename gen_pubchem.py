#!/usr/bin/env python


import os

import pubchempy
import pandas as pd
import click

import utils


def query_pubchem(formula, max_count=100):
    try:
        df = pubchempy.get_compounds(
                formula,
                "formula",
                listkey_count=max_count,
                as_dataframe=True
                )
        return df
    except:
        print(f"{formula} failed to query!")


def exhaust_query(formula, max_count=100):
    """
    This function will persist in querying until no more
    molecules are known in the database.
    """
    batch = list()
    pubchem_df = query_pubchem(formula, max_count)
    if pubchem_df is not None:
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
    """
    Function to generate XYZ files from a Pandas
    DataFrame from pubchem.

    Because PubChem also returns isotopologues,
    we will also prefilter the dataframe to remove
    them, and use the isomeric smiles instead of
    canonical smiles.
    """
    # Write smiles to file
    subset = pubchem_df.loc[
        ~pubchem_df["iupac_name"].str.contains("deut")
        ]
    # Get only unique smiles entries, so we don't do
    # redundant calculations
    smi = subset["isomeric_smiles"].unique()
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


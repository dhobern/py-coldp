from coldp import COLDP
import sys
import csv
import re
import pandas as pd
import numpy as np

pd.options.mode.chained_assignment = None  # default='warn'


def add_rank(names, rank):
    names[rank] = np.nan
    names.loc[names["rank"] == rank, rank] = names["scientificName"]
    loop = True
    while loop:
        names = pd.merge(
            names,
            names[~names[rank].isnull()][["ID_t", rank]],
            left_on="parentID",
            right_on="ID_t",
            how="left",
            suffixes=("", "_p"),
        )
        parent_rank = f"{rank}_p"
        count = len(names[(names[rank].isnull()) & (~names[parent_rank].isnull())])
        if count > 0:
            names.loc[
                (names[rank].isnull()) & (~names[parent_rank].isnull()), rank
            ] = names[parent_rank]
        else:
            loop = False
        names.drop(columns=["ID_t_p", parent_rank], inplace=True)
    return names


def build_hierarchy(coldp):
    names = pd.merge(
        coldp.names[["ID", "scientificName", "rank"]],
        coldp.taxa[["nameID", "ID", "parentID"]],
        left_on="ID",
        right_on="nameID",
        how="left",
        suffixes=("", "_t"),
    )
    accepted = names[~names["ID_t"].isnull()]
    accepted["status"] = "accepted"
    synonyms = names[names["ID_t"].isnull()]
    synonyms["status"] = "synonym"
    hierarchy_ranks = [
        "order",
        "superfamily",
        "family",
        "subfamily",
        "tribe",
        "subtribe",
        "genus",
        "subgenus",
        "species",
        "subspecies",
    ]
    for rank in hierarchy_ranks:
        accepted = add_rank(accepted, rank)
    synonyms = pd.merge(
        synonyms,
        coldp.synonyms[["nameID", "taxonID"]],
        left_on="ID",
        right_on="nameID",
        how="left",
        suffixes=("", "_s"),
    )
    synonyms = pd.merge(
        synonyms,
        accepted,
        left_on="taxonID",
        right_on="ID_t",
        how="left",
        suffixes=("", "_t"),
    )
    accepted.drop(columns=["nameID"], inplace=True)
    accepted.rename(columns={"ID": "nameID", "ID_t": "taxonID"}, inplace=True)
    synonyms.drop(
        columns=[
            "nameID",
            "parentID",
            "nameID_s",
            "taxonID",
            "ID_t",
            "scientificName_t",
            "rank_t",
            "nameID_t",
            "status_t",
        ],
        inplace=True,
    )
    synonyms.rename(
        columns={"ID": "nameID", "ID_t_t": "taxonID", "parentID_t": "parentID"},
        inplace=True,
    )
    hierarchy = pd.concat([accepted, synonyms], axis=0)
    hierarchy = hierarchy.fillna("ZZZ")
    hierarchy.loc[hierarchy["scientificName"] == "Lepidoptera", "superfamily"]
    hierarchy["sortkey"] = (
        hierarchy["order"]
        + hierarchy["superfamily"]
        + hierarchy["family"]
        + hierarchy["subfamily"]
        + hierarchy["tribe"]
        + hierarchy["subtribe"]
        + hierarchy["genus"]
        + hierarchy["subgenus"]
        + hierarchy["species"]
        + hierarchy["subspecies"]
    )
    hierarchy["sortkey"] = (
        hierarchy["sortkey"].apply(lambda x: re.sub(r"Z*$", "", str(x)))
        + hierarchy["status"].apply(lambda x: "AAA" if x == "accepted" else "ZZZ")
        + hierarchy["scientificName"]
    )
    hierarchy.sort_values("sortkey", inplace=True)
    hierarchy.drop(columns=["sortkey"], inplace=True)
    hierarchy.replace("ZZZ", "", inplace=True)
    return hierarchy


if __name__ == "__main__":
    if len(sys.argv) > 2:
        with open("show.csv", "w", encoding="utf8", newline="") as f:
            writer = csv.writer(f, delimiter=",")

            coldp = COLDP(sys.argv[1], sys.argv[2], issues_to_stdout=True)
            hierarchy = build_hierarchy(coldp)
            hierarchy = pd.merge(
                hierarchy.drop(columns=["scientificName", "rank"]),
                coldp.names.drop(
                    columns=["genus", "status", "remarks", "code", "link"]
                ),
                left_on="nameID",
                right_on="ID",
                how="left",
            ).drop("ID")
            hierarchy.to_csv("hierarchy.csv", index=False)
    else:
        print(f"Usage: python {sys.argv[0]} folder coldp_name")

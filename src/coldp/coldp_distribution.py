import os
import sys
import pandas
import numpy
import re
import shutil
import pathlib
from datetime import datetime
import warnings

#warnings.filterwarnings("ignore", category=numpy.VisibleDeprecationWarning)


# Return taxon id for taxon record or for taxon record linked via synonym record to a name record based on search_name
#   - supplied_name: the original name supplied by the caller
#   - search_name: original name or gender-flipped variant for actual search on this invocation
#   - taxa, names, synonyms: data frames from COLDP package
def get_taxon_id(supplied_name, search_name, taxa, names, synonyms):
    taxon_id = None
    taxon_name = None

    taxon = taxa[taxa["species"] == search_name]
    if len(taxon) == 1:
        taxon_id = taxon.iloc[0]["ID"]
        taxon_name = names[names["ID"] == taxon.iloc[0]["nameID"]].iloc[0][
            "scientificName"
        ]
    else:
        # Search name may be considered a synonym within the COLDP package
        name = names[names["scientificName"] == search_name]
        if len(name) == 1:
            synonym = synonyms[synonyms["nameID"] == name.iloc[0]["ID"]]
            if len(synonym) == 1:
                taxon = taxa[taxa["ID"] == synonym.iloc[0]["taxonID"]]
                if len(taxon) == 1:
                    taxon_id = taxon.iloc[0]["ID"]
                    taxon_name = names[names["ID"] == taxon.iloc[0]["nameID"]].iloc[0][
                        "scientificName"
                    ]
            else:
                # One remaining possibility is that there was more than 1 matching taxon because one or more
                # was for a subspecies - this solves these.
                taxon = taxa[taxa["nameID"] == name.iloc[0]["ID"]]
                if len(taxon) == 1:
                    taxon_id = taxon.iloc[0]["ID"]
                    taxon_name = names[names["ID"] == taxon.iloc[0]["nameID"]].iloc[0][
                        "scientificName"
                    ]

    # If the taxon_id has not been found and we are not already testing a gender-flipped name, try one
    if taxon_id is None and supplied_name == search_name:
        if supplied_name.endswith("us"):
            search_name, taxon_id, taxon_name = get_taxon_id(
                supplied_name,
                supplied_name[0 : len(supplied_name) - 2] + "a",
                taxa,
                names,
                synonyms,
            )
        elif supplied_name.endswith("a"):
            search_name, taxon_id, taxon_name = get_taxon_id(
                supplied_name,
                supplied_name[0 : len(supplied_name) - 1] + "us",
                taxa,
                names,
                synonyms,
            )

    return search_name, taxon_id, taxon_name


def validate_listitem(species_list, index, taxa, names, synonyms):
    taxon_id = None
    supplied_name = species_list.iloc[index]["SuppliedName"]
    search_name, taxon_id, taxon_name = get_taxon_id(
        supplied_name, supplied_name, taxa, names, synonyms
    )
    species_list.loc[
        species_list["SuppliedName"] == supplied_name,
        ["SearchName", "TaxonID", "TaxonName"],
    ] = [search_name, taxon_id, taxon_name]
    if "RegionList" in species_list.columns:
        return set(species_list.iloc[index]["RegionList"].split("|"))
    return set()


region_id = None
species_filename = None
coldp_foldername = None
data_foldername = None
reference_id = None
show_help = False

numeric = re.compile("^[0-9]+$")
iso2 = re.compile("^[A-Z][A-Z](-[A-Z0-9]+)?$")

# Be tolerant of arguments in any order - one should be a country code, one should be a reference ID, one should be a
# COLDP folder and one should be a CSV file
for i in range(1, len(sys.argv)):
    arg = sys.argv[i]
    if os.path.exists(arg):
        if os.path.isfile(arg) and species_filename is None:
            species_filename = arg
            print(f"Distribution list: {species_filename}")
        elif os.path.isdir(arg) and coldp_foldername is None:
            coldp_foldername = arg
            if os.path.isdir(os.path.join(coldp_foldername, "data")):
                data_foldername = os.path.join(arg, "data")
            else:
                data_foldername = coldp_foldername
            print(f"COLDP folder: {data_foldername}")
        else:
            show_help = True
    elif numeric.search(arg):
        reference_id = arg
        print(f"Reference id: {reference_id}")
    elif iso2.search(arg):
        region_id = arg
        print(f"Country: {region_id}")
    else:
        show_help = True

if show_help:
    print(
        "Usage: python "
        + sys.argv[0]
        + " <coldp_foldername> <species list file> <ISO 2-letter country code> <reference id>"
    )
    exit()

region_ids = set()
if region_id != "XX":
    region_ids.add(region_id)

print("Loading distributions")
distribution_filename = os.path.join(data_foldername, "distribution.csv")
distributions = None
if os.path.exists(distribution_filename):
    print("Loading distributions")
    distributions = pandas.read_csv(
        distribution_filename, dtype=str, keep_default_na=False
    )
if distributions is None:
    print("Could not get distributions")
    exit()

reference_filename = os.path.join(data_foldername, "reference.csv")
reference = None
if os.path.exists(reference_filename):
    print("Loading references")
    references = pandas.read_csv(reference_filename, dtype=str, keep_default_na=False)
    references = references[references["ID"] == reference_id]
    if len(references) == 1:
        reference = references.iloc[0]
if reference is None:
    print("Could not get reference")
    exit()

taxon_filename = os.path.join(data_foldername, "taxon.csv")
taxa = None
if os.path.exists(taxon_filename):
    print("Loading taxa")
    taxa = pandas.read_csv(taxon_filename, dtype=str, keep_default_na=False)
if taxa is None:
    print("Could not get taxa")
    exit()

name_filename = os.path.join(data_foldername, "name.csv")
names = None
if os.path.exists(name_filename):
    print("Loading names")
    names = pandas.read_csv(name_filename, dtype=str, keep_default_na=False)
if names is None:
    print("Could not get names")
    exit()

synonym_filename = os.path.join(data_foldername, "synonym.csv")
synonyms = None
if os.path.exists(synonym_filename):
    print("Loading synonyms")
    synonyms = pandas.read_csv(synonym_filename, dtype=str, keep_default_na=False)
if synonyms is None:
    print("Could not get synonyms")
    exit()

print("Loading species list")
species_list = pandas.read_csv(
    species_filename, header=None, dtype=str, keep_default_na=False
)

if region_id == "XX" and len(species_list.columns) < 3:
    print("Species list does not include regions")
    exit()

list_headings = ["SuppliedName", "Remarks", "RegionList"]
species_list.columns = list_headings[0 : len(species_list.columns)]
species_list["SearchName"] = species_list.iloc[:, 0]
species_list["TaxonID"] = numpy.nan
species_list["TaxonName"] = numpy.nan

for s in range(len(species_list)):
    region_ids.update(validate_listitem(species_list, s, taxa, names, synonyms))

missing_taxa = species_list[species_list["TaxonID"].isnull()]
if len(missing_taxa) > 0:
    print("\nFix the following missing taxa before proceeding:\n")
    for index, row in missing_taxa.iterrows():
        print("    " + row["SuppliedName"])
    exit()

region_filename = os.path.join(data_foldername, "region.csv")
bad_regions = []
if os.path.exists(region_filename):
    regions = pandas.read_csv(region_filename, dtype=str, keep_default_na=False)
    if regions is not None:
        for r in region_ids:
            matches = regions[regions["ID"] == r]
            if len(matches) != 1:
                bad_regions.append(r)

if len(bad_regions) > 0:
    print(
        "Could not find region"
        + (
            ("s:" + ", ".join(bad_regions))
            if len(bad_regions) != 1
            else (": " + bad_regions[0])
        )
    )
    exit()

print(
    "\nReference: "
    + reference["author"]
    + ", "
    + reference["issued"]
    + ", "
    + reference["title"]
)
print("\nRegions: ")
for r in region_ids:
    region = regions[regions["ID"] == r].iloc[0]
    print("    " + region["ID"] + " (" + region["name"] + ")")
print("\nSpecies: ")
for index, row in species_list.iterrows():
    if row["SuppliedName"] == row["TaxonName"]:
        print("    " + row["SuppliedName"] + " -> " + row["TaxonID"])
    else:
        print(
            "    "
            + row["SuppliedName"]
            + " -> "
            + row["SearchName"]
            + " -> "
            + row["TaxonName"]
            + " -> "
            + row["TaxonID"]
        )
        if "Remarks" in row and row["Remarks"]:
            row["Remarks"] = row["Remarks"] + " (as " + row["SuppliedName"] + ")"
        else:
            row["Remarks"] = "As " + row["TaxonName"]

print("\nProceed (Y/N)?")
response = input()
if response.lower() != "y":
    print("\nExiting")
    exit()

shutil.make_archive(
    os.path.join(
        pathlib.Path(coldp_foldername).parent.absolute(),
        "data-" + datetime.now().strftime("%Y-%m-%dT%H.%M.%S.%f"),
    ),
    "zip",
    data_foldername,
)

for r in region_ids:
    for index, item in species_list.iterrows():
        if region_id != "XX" or r in item["RegionList"]:
            taxon_id = item["TaxonID"]
            remark = item["Remarks"] if "Remarks" in item else ""
            existing = distributions[
                (distributions["taxonId"] == taxon_id)
                & (distributions["area"] == r)
                & (distributions["referenceID"] == "")
            ]
            if len(existing) == 1:
                distributions.loc[
                    (distributions["taxonId"] == taxon_id)
                    & (distributions["area"] == r)
                    & (distributions["referenceID"] == ""),
                    ["referenceID", "remarks"],
                ] = [reference_id, remark]
            else:
                existing = distributions[
                    (distributions["taxonId"] == taxon_id)
                    & (distributions["area"] == r)
                    & (distributions["referenceID"] == reference_id)
                ]
                if len(existing) == 0:
                    distributions.loc[len(distributions)] = [
                        taxon_id,
                        r,
                        "iso",
                        "native",
                        reference_id,
                        remark,
                    ]

distributions.to_csv(distribution_filename, index=False)

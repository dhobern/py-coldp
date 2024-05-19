#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Donald Hobern, dhobern@gmail.com
# Created Date: 2024-05-18
# version ='2024.0.2'
# ---------------------------------------------------------------------------
""" 
coldp_add is a Python tool for making batch updates to a Catalogue of Life Data 
Package (COLDP) archive using the :py:class:`~coldp.COLDP` class.

The initial version of this tool provides a simplified mechanism for adding 
wholly new species to a COLDP archive. Other behaviours (new genera and higher
ranks, new or modified synonymy and classification, additions to distribution
data, etc.) will be progressively added.

Usage: python coldp_add.py COLDP_path addition_filepath

* COLDP_path is the path to the root folder of a COLDP archive - the final part of the path is the name of the archive (e.g. if COLDP_path is "./archives/mammals", there should be a "mammals" folder in the folder "./archives", and the mammals folder should contain CSV/TSV files for the COLDP record types, optionally nested in a data subfolder)
* addition_filepath is the path to a CSV file with column names matching property names from COLDP data classes. Each record is a denormalised view of the data for a name and taxon. Where applicable, column names may include a lowercase namespace prefix indicating wthe record type to which the property belongs.

coldp_app will infer which elements relate to a reference and will find an existing record for this reference in the COLDP archive or will create a new reference. In either case, the ID for the reference will be included as the referenceID in all other COLDP records based on this row.

coldp_app will then use name elements to create a new species name record and taxon record in the appropriate genus inside the COLDP archive. IDs created for the name and taxon records will be used as nameID and taxonID as appropriate in other COLDP records based on this row.

coldp_app will then create any required typematerial, speciesinteraction or distribution records to accommodate remaining columns in the source record.
"""
from coldp import COLDP
import os
import sys
import shutil
from datetime import datetime
import pandas as pd

"""
Process command-line arguments
 
Folder + Package
Addition file
"""

argv_valid = True

if len(sys.argv) != 3:
    argv_valid = False
elif not os.path.isdir(sys.argv[1]):
    print(f"COLDP folder not found {sys.argv[1]}\n")
    argv_valid = False
elif not os.path.isfile(sys.argv[2]):
    print(f"COLDP addition file not found {sys.argv[2]}\n")
    argv_valid = False

if argv_valid:
    df = pd.read_csv(sys.argv[2], dtype=str, keep_default_na=False)
    if df is None:
        print(f"COLDP addition file not readable as CSV {sys.argv[2]}\n")
        argv_valid = False

if argv_valid:
    path, package = os.path.split(sys.argv[1])
    coldp = COLDP(package, path, issues_to_stdout=True)
    if coldp is None:
        print(f"Folder not readable as COLDP package {sys.argv[1]}\n")
        argv_valid = False

if not argv_valid:
    print(f"Usage: python {sys.argv[0]} COLDP_path addition_filepath")
    exit(0)

backup_file = f"{sys.argv[1]}-{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"
print(f"{sys.argv[1]}: Archiving to {backup_file}")
shutil.make_archive(backup_file, "zip", sys.argv[1])


headings = coldp.get_available_column_headings()
warned = []

for index, row in df.iterrows():
    records = {}
    for c in df.columns:
        if ":" in c:
            nsc = c
            t, c = c.split(":")[0:2]
            if t in headings and c in headings[t]:
                if t in records:
                    r = records[t]
                else:
                    r = {}
                    records[t] = r
                r[c] = row[nsc]
            else:
                print(f"Could not handle column {t} : {c}")
                exit()
        else:
            found = 0
            for t in headings.keys():
                if c in headings[t]:
                    if found > 0 and c not in warned:
                        print(
                            f"Warning: column {c} found to multiple COLDP tables - supplied values will appear in all matching tables"
                        )
                        warned.append(c)
                    found += 1
                    if t in records:
                        r = records[t]
                    else:
                        r = {}
                        records[t] = r
                    r[c] = row[c]
            if found == 0:
                print(f"Could not handle column {c}")
                exit(0)

    if "reference" in records:
        records["reference"] = coldp.add_reference(records["reference"])
        print(f"Reference: {records['reference']}")
        for t in headings:
            if t in records and "referenceID" in headings[t]:
                records[t]["referenceID"] = records["reference"]["ID"]

    taxonID = None

    if "name" in records:
        name = records["name"]
        bundle = coldp.start_name_bundle(name)
        parentID = None
        if "rank" in name and name["rank"] == "species":
            if "genus" in name:
                print(f"Looking for genus {name['genus']}")
                parent = coldp.find_taxon(name["genus"], None, "genus")
                if parent is not None:
                    parentID = parent["ID"]
                    print(f"Parent ID: {parentID}")
            else:
                print(f"Cannot find genus element in name record {name}")
        coldp.add_names(bundle, parentID)
        print(f"Name: {bundle.accepted}")
        print(f"Taxon ID: {bundle.accepted_taxon_id}")
        for t in headings:
            if t in records:
                if "nameID" in headings[t]:
                    records[t]["nameID"] = bundle.accepted["ID"]
                if "taxonID" in headings[t]:
                    taxonID = bundle.accepted_taxon_id
                    records[t]["taxonID"] = taxonID

    if "taxon" in records and taxonID is not None:
        coldp.modify_taxon(taxonID, records["taxon"])

    if "typematerial" in records:
        print(f"Type material: {records['typematerial']}")
        coldp.add_type_material(records["typematerial"])

    if "speciesinteraction" in records:
        print(f"Species interaction: {records['speciesinteraction']}")
        coldp.add_species_interaction(records["speciesinteraction"])

    if "distribution" in records:
        print(f"Distribution: {records['distribution']}")
        coldp.add_distribution(records["distribution"])

coldp.save()

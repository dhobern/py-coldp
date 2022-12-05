import os
import sys
import csv
import pandas as pd
import numpy as np
import re
import shutil
import pathlib
from datetime import datetime

source_foldername = None
target_foldername = None
force = False

name_from_nameusage = { 
    "status": "nameStatus",
    "genus": "genericName",
    "referenceID": "nameReferenceID",
    "remarks": "nameRemarks"
}

taxon_from_nameusage = { 
    "nameID": "ID",
    "provisional": "status"
}

synonym_from_nameusage = {
    "taxonID": "parentID",
    "nameID": "ID"
}

nameusage_columns = [ 
    "ID", "sourceID", "parentID", "basionymID", "status", "scientificName", 
    "authorship", "rank", "notho", "uninomial", "genericName", 
    "infragenericEpithet", "specificEpithet", "infraspecificEpithet", 
    "cultivarEpithet", "namePhrase", "nameReferenceID", "publishedInYear", 
    "publishedInPage", "publishedInPageLink", "code", "nameStatus", 
    "accordingToID", "accordingToPage", "accordingToPageLink", "referenceID", 
    "scrutinizer", "scrutinizerID", "scrutinizerDate", "extinct", 
    "temporalRangeStart", "temporalRangeEnd", "environment", "species", 
    "section", "subgenus", "genus", "subtribe", "tribe", "subfamily", "family", 
    "superfamily", "suborder", "order", "subclass", "class", "subphylum", 
    "phylum", "kingdom", "sequenceIndex", "branchLength", "link", "nameRemarks", 
    "remarks" ]

name_columns = [ 
    "ID", "sourceID", "basionymID", "scientificName", 
    "authorship", "rank", "uninomial", "genus", 
    "infragenericEpithet", "specificEpithet", "infraspecificEpithet", 
    "cultivarEpithet", "code", "status", "referenceID", "publishedInYear", 
    "publishedInPage", "publishedInPageLink", "link", "remarks" ]

taxon_columns = [
    "ID", "sourceID", "parentID", "nameID", "namePhrase", "nameReferenceID", 
    "publishedInYear", "publishedInPage", "publishedInPageLink", "code", 
    "nameStatus", "accordingToID", "scrutinizer", "scrutinizerID", 
    "scrutinizerDate", "referenceID", "extinct", "temporalRangeStart", 
    "temporalRangeEnd", "environment", "species", "section", "subgenus", 
    "genus", "subtribe", "tribe", "subfamily", "family", 
    "superfamily", "suborder", "order", "subclass", "class", "subphylum", 
    "phylum", "kingdom", "link", "remarks" ]

rank_columns = [ "species", "section", "subgenus", "genus", "subtribe", "tribe", 
    "subfamily", "family", "superfamily", "suborder", "order", "subclass", 
    "class", "subphylum", "phylum", "kingdom" ]

synonym_columns = [
    "ID", "sourceID", "taxonID", "nameID", "namePhrase", "accordingToID", 
    "status", "referenceID",  "link", "remarks" ]

def read_table(folder, name, extension):
    if extension in [ "tsv" ]:
        delimiter = "\t"
    elif extension in [ "csv", "txt"]:
        delimiter = ","
    else:
        print(f"Unexpected extension: {name}.{extension}")
        delimiter = None
        
    if delimiter is not None:
        table = pd.read_csv(os.path.join(folder, name), dtype=str, 
                keep_default_na=False, sep='\t')
        columns = table.columns.tolist()
        for i in range(len(columns)):
            if ":" in columns[i]:
                columns[i] = columns[i].split(":")[1]
        table.columns = columns
    else:
        table = None
    
    return table

def open_writer(folder, name, headings):
    file = open(os.path.join(folder, name), 'w', newline='', encoding="utf8")
    writer = csv.writer(file, delimiter=',')
    writer.writerow(headings)
    return writer

def fix_basionyms(names, synonyms):
    for index, row in names.iterrows():
        if row["authorship"].startswith("("):
            authorship = row["authorship"][1:len(row["authorship"])-1]
            if row["infraspecificEpithet"]:
                epithet = row["infraspecificEpithet"]
            else:
                epithet = row["specificEpithet"]
            match = names[(names["authorship"] == authorship) & 
                    ((names["infraspecificEpithet"] == epithet) | 
                        ((names["infraspecificEpithet"] == "") & 
                            (names["specificEpithet"] == epithet)))]
            if len(match) > 1:
                synonym_match = pd.merge(synonyms[synonyms["taxonID"] == row["ID"]][["nameID", "taxonID"]], names, left_on="nameID", right_on="ID", how="left")[["ID"]]
                match = pd.merge(match, synonym_match, on="ID", how="inner")
            if len(match) == 0:
                print(f"Cannot find basionym for {row['ID']}: {row['scientificName']} {row['authorship']}")
            elif len(match) == 1:
                names.loc[names["ID"] == row["ID"], "basionymID"] = match.iloc[0]["ID"]
        else:
            names.loc[names["ID"] == row["ID"], "basionymID"] = row["ID"]

def fix_classification(taxa, ranks, parent):
    if parent is None:
        parents = pd.merge(taxa, ranks[ranks["rank"] == "kingdom"], 
                on="nameID", how="right")
    else:
        '''
        classification = []
        for i in range(len(rank_columns)):
            classification.append(parent[rank_columns[i]])
        '''
        classification = parent[rank_columns]
        if parent["rank"] in rank_columns:
            classification[rank_columns.index(parent["rank"])] = parent["scientificName"]

        taxa.loc[taxa["parentID"] == parent["ID"], rank_columns] \
                = classification.to_numpy()
        #print(",".join(taxa.loc[taxa["parentID"] == parent["ID"], [rank_columns]].array))
        parents = pd.merge(taxa.loc[taxa["parentID"] == parent["ID"]],
                ranks, on="nameID", how="left")
    for i in range(len(parents)):
        fix_classification(taxa, ranks, parents.iloc[i])

def write_row(table, row, columns, mappings):
    out_row = {}
    for column in columns:
        #print(f"Column: {column}")
        if column in mappings.keys():
            mapped = mappings[column]
            #print(f"Mapped: {mapped}")
        else:
            mapped = column
        if mapped in row and row[mapped]:
            #print(f"Value: {row[mapped]}")
            out_row[column] = row[mapped]
        else:
            out_row[column] = ""
    return table.append(out_row, ignore_index=True)

def save_table(table, folder, name):
    if len(table) > 0:
        table.replace('', np.nan, inplace=True)
        table.dropna(how='all', axis=1, inplace=True)
        table.replace(np.nan, '', inplace=True)
        table.to_csv(os.path.join(folder, name + ".csv"), index=False)

for i in range(1, len(sys.argv)):
    if sys.argv[i].startswith("-"):
        if sys.argv[i][1:].lower() in [ "force", "f"]:
            force = True
        else:
            target_foldername = None
            break
    elif not source_foldername:
        if os.path.isdir(sys.argv[i]):
            source_foldername = sys.argv[i]
        else:
            break
    elif not target_foldername and sys.argv[i] != source_foldername:
        target_foldername = sys.argv[i]
    else:
        target_foldername = None
        break

if not target_foldername:
    print(f"python {sys.argv[0]} [-F[orce]] source_foldername target_foldername")
    sys.exit()

target_basefoldername = target_foldername
if os.path.basename(target_foldername) != "data":
    target_foldername = os.path.join(target_foldername, "data")

if not os.path.exists(target_foldername):
    if force:
        os.makedirs(target_foldername)
    else:
        print(f"target_foldername must exist or -Force must be specified")
elif force:
    for filename in os.listdir(target_foldername):
        file_path = os.path.join(target_foldername, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

pd.set_option('display.max_columns', None)

for filename in os.listdir(source_foldername):
    parts = filename.lower().split(".")
    if len(parts) == 2:
        table_name = parts[0]
        if table_name == "nameusage":
            nameusages_in = read_table(source_foldername, filename, parts[1])
            names_out = pd.DataFrame(columns=name_columns)
            taxa_out = pd.DataFrame(columns=taxon_columns)
            synonyms_out = pd.DataFrame(columns=synonym_columns)
            for index, row in nameusages_in.iterrows():
                names_out = write_row(names_out, row, name_columns, 
                        name_from_nameusage)
                if row["status"] in \
                        [ "synonym", "ambiguous synonym", "misapplied" ]:
                    synonyms_out = write_row(synonyms_out, row, synonym_columns, 
                            synonym_from_nameusage)
                else:
                    taxa_out = write_row(taxa_out, row, taxon_columns, 
                            taxon_from_nameusage)
            ranks = names_out.loc[:, ["ID", "rank", "scientificName"]]
            ranks.columns = [ "nameID", "rank", "scientificName" ]
            fix_classification(taxa_out, ranks, None)
            fix_basionyms(names_out, synonyms_out)
            save_table(names_out, target_foldername, "name")
            save_table(taxa_out, target_foldername, "taxon")
            save_table(synonyms_out, target_foldername, "synonym")
        elif parts[1] in [ "tsv", "csv", "txt" ]:
            table_in = read_table(source_foldername, filename, parts[1])
            save_table(table_in,target_foldername, parts[0])
        else:
            shutil.copy(os.path.join(source_foldername, filename), 
                    target_basefoldername)
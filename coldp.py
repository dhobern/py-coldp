#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  : Donald Hobern, dhobern@gmail.com
# Created Date: 2022-11-27
# version ='1.0'
# ---------------------------------------------------------------------------
""" 
COLDP.py

COLDP is a Python class for creating or loading, manipulating and serialising 
Catalogue of Life Data Package (COLDP) files.

COLDP is the preferred format within Catalogue of Life for organising taxonomic 
checklist datasets. See: https://github.com/CatalogueOfLife/coldp 
"""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import pandas as pd
import numpy as np  
import logging
import os
import re
import shutil

# ---------------------------------------------------------------------------
# Column headings for each data class within COLDP packages
# ---------------------------------------------------------------------------
nameusage_headings = [ 
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

name_headings = [ 
    "ID", "sourceID", "basionymID", "scientificName", 
    "authorship", "rank", "uninomial", "genus", 
    "infragenericEpithet", "specificEpithet", "infraspecificEpithet", 
    "cultivarEpithet", "code", "status", "referenceID", "publishedInYear", 
    "publishedInPage", "publishedInPageLink", "link", "remarks" ]

taxon_headings = [
    "ID", "sourceID", "parentID", "nameID", "namePhrase", "nameReferenceID", 
    "publishedInYear", "publishedInPage", "publishedInPageLink", "code", 
    "nameStatus", "accordingToID", "scrutinizer", "scrutinizerID", 
    "scrutinizerDate", "referenceID", "extinct", "temporalRangeStart", 
    "temporalRangeEnd", "environment", "species", "section", "subgenus", 
    "genus", "subtribe", "tribe", "subfamily", "family", 
    "superfamily", "suborder", "order", "subclass", "class", "subphylum", 
    "phylum", "kingdom", "link", "remarks" ]

rank_headings = [ "kingdom", "phylum", "subphylum", "class", "subclass", 
    "order", "suborder", "superfamily", "family", "subfamily", "tribe", 
    "subtribe", "genus", "subgenus", "section", "species" ]

synonym_headings = [
    "ID", "sourceID", "taxonID", "nameID", "namePhrase", "accordingToID", 
    "status", "referenceID",  "link", "remarks" ]

reference_headings = [ 
    "ID", "author", "title", "issued", "containerTitle", "volume", "issue", 
    "page", "link", "citation" ]

type_material_headings = [ 
    "nameID", "citation", "status", "referenceID", "locality", "country", 
    "latitude", "longitude", "elevation", "date", "collector", 
    "institutionCode", "sex", "remarks" ]

name_relation_headings = [ 
    "nameID", "relatedNameID", "type", "referenceID", "remarks" ]

distribution_headings = [ 
    "taxonId", "area", "gazetteer", "status", "referenceID", "remarks" ]

species_interaction_headings = [ 
    "taxonID", "relatedTaxonID", "relatedTaxonScientificName", "type", 
    "referenceID", "remarks" ]

# ---------------------------------------------------------------------------
# Column headings for the issues log
# ---------------------------------------------------------------------------
issues_headings = [ "context", "issue" ]

# ---------------------------------------------------------------------------
# Dictionaries for mapping columns from nameusage format into separate
# name, taxon and synonym dataframes. 
# 
# Keys indicate the column name in the target table and values represent the
# column name in the nameusage table 
# ---------------------------------------------------------------------------
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

#----------------------------------------------------------------------------
# Dictionary for mappings between id columns in tables
#
# The keys to the outer dictionary are table names (for use with 
# table_by_name()) which contain an "ID" column. Each key is associated with
# an inner dictionary that maps the name of another (or sometimes the same
# table) to a list of columns that contain ID values from the first table. 
#----------------------------------------------------------------------------
id_mappings = { "name": { "name": ["basionymID"],
                       "taxon": ["nameID"],
                       "name_relation": ["nameID", "relatedNameID"],
                       "synonym": ["nameID"] },
                "taxon": { "taxon": ["parentID"],
                       "distribution": ["taxonID"],
                       "species_interaction": ["taxonID", "relatedTaxonID"],
                       "synonym": ["taxonID"] },
                "reference": { "name": ["sourceID", "referenceID"],
                       "taxon": ["sourceID", "referenceID", "nameReferenceID"],
                       "synonym": ["sourceID", "referenceID"],
                       "type_material": ["referenceID"],
                       "distribution": ["referenceID"],
                       "species_interaction": ["referenceID"],
                       "name_relation": ["referenceID"] },
                "synonym": {} }

# ---------------------------------------------------------------------------
# Recognised file extensions and their associated separator characters
# ---------------------------------------------------------------------------
csv_extensions = { 
    "csv" : ",", 
    "tsv": "\t", 
    "txt": "\t" 
}

# ---------------------------------------------------------------------------
# Regular expression patterns for uninomial and species-rank epithets and a
# superset expression including both
# ---------------------------------------------------------------------------
uninomial_pattern = re.compile('^[A-Z][a-z]+$')
epithet_pattern = re.compile('^[a-z]-?[a-z]+$')
name_pattern = re.compile('^[A-Za-z]-?[a-z]+$')

# ---------------------------------------------------------------------------
# 
# CLASS: NameBundle
#
# Wrapper class to manage set of associated names, normally an accepted name
# and one or more synonyms. The bundle handles the logic of associated name 
# variations.
# 
# Usage:
#   The minimal usage is to create a new bundle with an accepted name. One 
#   or more synonyms can also be supplied using the add_synonym() method.
#   The bundle is then submitted to the COLDP add_names method for 
#   processing.
#
# ---------------------------------------------------------------------------
class NameBundle:

    #------------------------------------------------------------------------
    # __init__ - Initialise NameBundle
    #
    #   coldp: COLDP object to handle logging of any issues and normalise 
    #       name records
    #   accepted: Dictionary mapping COLDP name elements to values for the 
    #       accepted name - a name record and a taxon record will be added to 
    #       the COLDP package for the accepted name
    #   incertae_sedis: Flag to indicate if the resulting taxon record should
    #       be marked "incertae sedis" 
    #------------------------------------------------------------------------
    def __init__(self, coldp, accepted:dict, 
                incertae_sedis:bool=False, sic:bool=False) -> None:
        self.coldp = coldp
        if "rank" not in accepted:
            self.coldp.issue("Accepted name missing rank, assume species: " \
                    + str(accepted))
            accepted["rank"] = "species"
        self.accepted_sic = sic
        self.accepted = self.normalise_name(accepted, sic)
        self.incertae_sedis = incertae_sedis
        self.accepted_taxon_id = None
        self.synonyms = []
        self.species = None
        self.species_synonym = None
        self.accepted_species_id = None

    #------------------------------------------------------------------------
    # add_synonym - Register an additional name as a synonym for the accepted
    #       name
    #
    #   synonym: Dictionary mapping COLDP name elements to values for the 
    #       synonymous name - a name record and a synonym record will be 
    #       added to the COLDP package for the synonym
    #   sic: Flag to indicate that the name does not follow expected 
    #       formatting rules for a code-compliant name and that no issues
    #       should be logged for this 
    #------------------------------------------------------------------------
    def add_synonym(self, synonym:dict, sic:bool=False) -> None:
        normalised = self.normalise_name(synonym, sic)
        present = (normalised["scientificName"] == self.accepted["scientificName"])
        if not present:
            for s in self.synonyms:
                if normalised["scientificName"] == s["scientificName"]:
                    present = True
                    break
        if not present:
            self.synonyms.append(self.normalise_name(synonym, sic))

    #------------------------------------------------------------------------
    # normalise_name - Ensure that a name record dictionary contains all 
    #       necessary/appropriate elements
    #
    #   name: Dictionary containing name to be normalised
    #   sic: Flag to indicate that the name does not follow expected 
    #       formatting rules for a code-compliant name and that no issues
    #       should be logged for this 
    #
    #   Returns name dictionary updated with extra values 
    #------------------------------------------------------------------------
    def normalise_name(self, name:dict, sic:bool=False) -> dict:

        # Ensure main elements are present in dictionary
        for k in [
                "uninomial", "genus", "infragenericEpithet", 
                "specificEpithet", "infraspecificEpithet", "authorship", 
                "rank", "basionymID", "scientificName" ]:
            if k not in name:
                name[k] = ""

        # Skip processing if sic
        if not sic:

            # Check uninomials for formatting - log issues and fix case
            for k in [ "uninomial", "genus", "infragenericEpithet" ]:
                if len(name[k]) > 0:
                    if not name_pattern.match(name[k]):
                        self.coldp.issue(
                            "Invalid pattern for supraspecific name: \"" + \
                                name[k] + "\"")
                    elif name[k][0].islower():
                        self.coldp.issue("Lowercase initial letter " +\
                            "for supraspecific name: " + name[k])
                        name[k] = name[k][0].upper() + name[k][1:]

            # Check epithets for formatting - log issues and fix case
            for k in [ "specificEpithet", "infraspecificEpithet" ]:
                if len(name[k]) > 0:
                    if not name_pattern.match(name[k]):
                        self.coldp.issue(
                                "Invalid pattern for epithet: " + name[k])
                    elif name[k][0].isupper():
                        self.coldp.issue(
                            "Uppercase initial letter for epithet: " + \
                                name[k])
                        name[k] = name[k][0].lower() + name[k][1:]

        if not sic or not name["scientificName"]:
            # Generate properly formatted scientific name
            if self.coldp.is_species_group(name):
                name["scientificName"], name["rank"] = \
                    self.coldp.construct_species_rank_name(
                            name["genus"],
                            name["infragenericEpithet"],
                            name["specificEpithet"],
                            name["infraspecificEpithet"],
                            name["rank"])
            else:
                name["scientificName"] = name["uninomial"]
            
        # Return normalised version of name
        return name

    #------------------------------------------------------------------------
    # derive_name - Use supplied values to create new name dictionary with
    #       missing elements copied from an existing name dictionary
    #
    #   name: Dictionary containing name on which new name is to be based
    #   values: Dictionary of values to override values in name 
    #
    #   Returns name dictionary with supplied values supplemented from name
    #------------------------------------------------------------------------
    def derive_name(self, name:dict, values:dict, sic:bool=False) -> dict:
        if "authorship" not in values and "authorship" in name and \
                    not name["authorship"].startswith("("):
            if "genus" in values or "specificEpithet" in values:
                values["authorship"] = "(" + name["authorship"] + ")"

        for k in [ 
                "basionymID", "authorship", "rank", "uninomial", "genus", 
                "infragenericEpithet", "specificEpithet", 
                "infraspecificEpithet", "code", "status" ]:
            if k in name and k not in values:
                values[k] = name[k]
        return self.normalise_name(values, sic)

# ---------------------------------------------------------------------------
#
# CLASS: COLDP
#
#   Class to manage a set of Pandas dataframes for CSV tables in Catalogue 
#   of Life Data Package.
# 
# Usage:
#   Initialise the COLDP object with folder, COLDP package name and other 
#   options. If the package exists in the folder, the COLDP package is 
#   initialised as a set of dataframes based on its data. Otherwise empty
#   dataframes are created.
# 
#   Data is added using the add_references, add_names, add_typematerial,
#   add_distribution and modify_taxon methods. The start_name_bundle is
#   used to produce a NameBundle object for preparing an accepted name
#   and associated synonyms to pass to add_names.
# 
#   The save method writes the data back to the same or another folder. 
# 
# ---------------------------------------------------------------------------
class COLDP:

    def __init__(self, folder, name, **kwargs):
        self.default_taxon = {
                "provisional": "false",
                "extinct": "false",
                "temporalRangeEnd": "Holocene"
            }
        self.code = "ICZN"
        self.species_from_trinomials = False
        self.subspecies_from_trinomials = False
        self.subgenus_free_synonyms = False
        self.basionyms_from_synonyms = False
        self.classification_from_parents = False
        self.allow_repeated_binomials = False
        self.issues_to_stdout = False

        if kwargs:
            self.set_options(kwargs)
        self.folder = folder
        self.name = name
        if folder is not None:
            source_folder = os.path.join(folder, name)
        else:
            source_folder = None
        self.name_usages = None
        self.names = self.initialise_dataframe(
                source_folder, "name", name_headings)
        self.taxa = self.initialise_dataframe(
                source_folder, "taxon", taxon_headings)
        self.synonyms = self.initialise_dataframe(
                source_folder, "synonym", synonym_headings)
        self.references = self.initialise_dataframe(
                source_folder, "reference", reference_headings)
        self.distributions = self.initialise_dataframe(
                source_folder, "distribution", distribution_headings)
        self.species_interactions = self.initialise_dataframe(
                source_folder, "speciesinteration", 
                species_interaction_headings)
        self.name_relations = self.initialise_dataframe(
                source_folder, "namerelation", name_relation_headings)
        self.type_materials = self.initialise_dataframe(
                source_folder, "typematerial", type_material_headings)
        if self.classification_from_parents:
            self.fix_classification()
        if self.basionyms_from_synonyms:
            self.fix_basionyms(self.names, self.synonyms)
        self.issues = None
        self.context = None

    def set(self, **kwargs):
        self.set_options(kwargs)

    def set_options(self, options):
        for key, value in options.items():
            if key == "insert_species_for_trinomials":
                self.species_from_trinomials = value
            elif key == "create_subspecies_for_infrasubspecifics":
                self.subspecies_from_trinomials = value
            elif key == "create_synonyms_without_subgenus":
                self.subgenus_free_synonyms = value
            elif key == "default_taxon_record":
                self.default_taxon = value
            elif key == "context":
                self.context = value
            elif key == "issues_to_stdout":
                self.issues_to_stdout = value
            elif key == "basionyms_from_synonyms":
                self.basionyms_from_synonyms = value
            elif key == "classification_from_parents":
                self.classification_from_parents = value
            elif key == "allow_repeated_binomials":
                self.allow_repeated_binomials = value
            elif key == "code" and value in ["ICZN", "ICBN"]:
                self.code = value

    def set_context(self, context):
        self.context = context

    def initialise_dataframe(self, foldername, name, default_headings):
        df = None
        nameusage_filename = None
        if foldername is not None and os.path.isdir(foldername):
            if os.path.isdir(os.path.join(foldername, "data")):
                foldername = os.path.join(foldername, "data")
            for filename in os.listdir(foldername):
                parts = filename.lower().split(".")
                if parts[0] == "nameusage":
                    nameusage_filename = filename
                if len(parts) == 2 and parts[0] == name \
                    and parts[1] in csv_extensions.keys():
                    df = pd.read_csv(os.path.join(foldername, filename), 
                            dtype=str, keep_default_na=False, 
                            sep=csv_extensions[parts[1]])
                    columns = df.columns.tolist()
                    for i in range(len(columns)):
                        if ":" in columns[i]:
                            columns[i] = columns[i].split(":")[1]
                    df.columns = columns
                    break
        if df is None and nameusage_filename is not None and \
                name in [ "name", "taxon", "synonym" ]:
            if self.name_usages is None:
                self.name_usages = self.initialise_dataframe(
                        foldername, "nameusage", nameusage_headings)
            if self.name_usages is not None:
                synonym_statuses = [ 
                        "synonym", "ambiguous synonym", "misapplied" ]
                if name == "name":
                    df = self.extract_table(
                            self.name_usages, name_headings, name_from_nameusage)
                elif name == "synonym":
                    df = self.extract_table(
                        self.name_usages[
                                self.name_usages.status.isin(synonym_statuses)], 
                        synonym_headings, synonym_from_nameusage)
                else:
                    df = self.extract_table(self.name_usages[
                            ~self.name_usages.status.isin(synonym_statuses)], 
                        taxon_headings, taxon_from_nameusage)

        if df is None and default_headings is not None:
            df = pd.DataFrame(columns=default_headings)

        return df

    def extract_table(self, df, headings, mappings):
        table = None

        headings_out = []
        headings_in = []

        for h in headings:
            m = mappings[h] if h in mappings else h
            if m in df.columns:
                headings_out.append(h)
                headings_in.append(m)

        if len(headings_in) > 0:
            table = df.loc[:, headings_in]
            table.columns = headings_out

        return table

    def fix_basionyms(self, names, synonyms):
        for index, row in names.iterrows():
            if row["authorship"].startswith("("):
                authorship = self.get_original_authorship(row["authorship"])
                if row["infraspecificEpithet"]:
                    epithet = row["infraspecificEpithet"]
                else:
                    epithet = row["specificEpithet"]
                match = names[(names["authorship"] == authorship) & 
                        ((names["infraspecificEpithet"] == epithet) | 
                            ((names["infraspecificEpithet"] == "") & 
                                (names["specificEpithet"] == epithet)))]
                if len(match) > 1:
                    synonym_match = pd.merge( \
                        synonyms[synonyms["taxonID"] == row["ID"]] \
                                [["nameID", "taxonID"]], 
                        names, left_on="nameID", right_on="ID", 
                        how="left")[["ID"]]
                    match = pd.merge(match, synonym_match, on="ID", how="inner")
                if len(match) == 0:
                    print(f"Cannot find basionym for {row['ID']}: " + \
                            f"{row['scientificName']} {row['authorship']}")
                elif len(match) == 1:
                    names.loc[names["ID"] == row["ID"], "basionymID"] = \
                            match.iloc[0]["ID"]
            else:
                names.loc[names["ID"] == row["ID"], "basionymID"] = row["ID"]
        for index, row in names[(names["infraspecificEpithet"] != "") &
                            (names["specificEpithet"] 
                                == names["infraspecificEpithet"])].iterrows():
            species = names[(names["genus"] == row["genus"]) &
                                (names["authorship"] == row["authorship"]) &
                                (names["rank"] == "species")]
            if len(species) == 1:
                names.loc[names["ID"] == row["ID"], "basionymID"] = \
                        species.iloc[0]["basionymID"]

    def fix_classification(self):
        ranks = self.names.loc[:, ["ID", "rank", "scientificName"]]
        ranks.columns = [ "nameID", "rank", "scientificName" ]
        self.fix_classification_recursive(self.taxa, ranks, None)

    def fix_classification_recursive(self, taxa, ranks, parent):
        if parent is None:
            for p in range(len(rank_headings)):
                parents = pd.merge(
                        taxa, ranks[ranks["rank"] == rank_headings[p]], 
                        on="nameID", how="right")
                if parents is not None and len(parents) > 0:
                    break
        else:
            classification = parent[rank_headings]
            if parent["rank"] in rank_headings:
                classification[rank_headings.index(parent["rank"])] = \
                    parent["scientificName"]

            taxa.loc[taxa["parentID"] == parent["ID"], rank_headings] \
                    = classification.to_numpy()
            parents = pd.merge(taxa.loc[taxa["parentID"] == parent["ID"]],
                    ranks, on="nameID", how="left")
        for i in range(len(parents)):
            self.fix_classification_recursive(taxa, ranks, parents.iloc[i])

    def sort_taxa(self):
        df = pd.merge(self.taxa, self.names, left_on="nameID", right_on="ID", 
                how="left", suffixes=["", "_name"]) \
                    [["ID", "parentID", "rank", "scientificName"]] \
                        .sort_values(["scientificName"])

        ids = []
        for index, row in df[(df["parentID"] == "") | \
                    (df["parentID"].isnull())].iterrows():
            ids = self.sort_taxa_recursive(df, ids, row["ID"])
        ids = ids + df[~df["ID"].isin(np.array(ids))]["ID"].tolist()

        sort_index = dict(zip(ids, range(len(ids))))
        self.taxa["sequence"] = self.taxa["ID"].map(sort_index)
        self.taxa = self.taxa.sort_values(["sequence"])
        self.taxa.drop("sequence", axis=1, inplace=True)
        self.taxa = self.taxa.reset_index(drop=True)

    def sort_taxa_recursive(self, df, ids, id):
        #ids.append(str(len(ids) + 1))
        ids.append(id)
        for index, row in df[df["parentID"] == id].iterrows():
            ids = self.sort_taxa_recursive(df, ids, row["ID"])
        return ids

    def sort_names(self):
        taxon_subset = self.taxa[["ID", "nameID"]]
        taxon_subset["idx"] = self.taxa.index
        sort_table = pd.merge(self.names[["ID"]], taxon_subset, 
            left_on="ID", right_on="nameID", how="left", 
            suffixes=["", "_taxon"])
        sort_table = pd.merge(sort_table, self.synonyms[["taxonID", "nameID"]], 
            left_on="ID", right_on="nameID", how="left", 
            suffixes=["", "_synonym"])
        sort_table.loc[sort_table["idx"].isnull(), ["idx"]] = \
            pd.merge(sort_table.loc[sort_table["idx"].isnull()], 
            taxon_subset, left_on="taxonID", right_on="ID", how="left")["idx_y"]
        self.names[["idx"]] = sort_table[["idx"]]
        self.names = self.names.sort_values(["idx", "scientificName"])
        self.names.drop(["idx"], axis=1, inplace=True)
        self.names = self.names.reset_index(drop=True)

    def table_by_name(self, name):
        tables = {
            "name": self.names, "taxon": self.taxa, "synonym": self.synonyms,
            "reference": self.references, "type_material": self.type_materials,
            "distribution": self.distributions, 
            "species_interaction": self.species_interactions,
            "name_relation": self.name_relations }
        return tables[name] if name in tables.keys() else None

    def reset_ids(self, name=None, prefix=None):
        for table_name in [ "name", "taxon", "reference", "synonym" ]:
            if name is None or table_name == name:
                table = self.table_by_name(table_name)
                if prefix is None:
                    prefix = ""
                if prefix == "" and table_name == "synonym":
                    prefix = "s_"
                table["newID"] = prefix + (table.index + 1).astype(str)
                by_ids = table.set_index("ID")
                joins = id_mappings[table_name]
                for key in joins.keys():
                    other_table = self.table_by_name(key)
                    for value in joins[key]:
                        if value in other_table.columns:
                            other_table[value] = \
                                other_table[value].map(by_ids["newID"])
                table["ID"] = table["newID"]
                table.drop("newID", axis=1, inplace=True)

    #
    # Ensure one or more references are included in references dataframe
    # 
    # Parameters:
    #   reference_list - list of dictionaries - each dictionary contains 
    #       values keyed by terms from reference_headings 
    #
    # Behaviour:
    #   Find existing ID values for each supplied reference, based on identity
    #   of: author, title, issued, containerTitle, volume, issue and page. Add 
    #   ID to the appropriate reference dictionary in references. If none 
    #   found, set ID to next index and add to references  
    # 
    # Returns:
    #   Updated copy of reference_list 
    #
    def add_references(self, reference_list):
        for i in range(len(reference_list)):
            has_content = False
            for k in reference_list[i].keys():
                if reference_list[i][k]:
                    has_content = True
                    break
            if has_content:
                r = reference_list[i]
                reference = self.find_reference(reference_list[i])
                if reference is not None:
                    reference_list[i] = reference
                    logging.debug(
                        "Matched reference ID " + reference_list[i]["ID"])
                else:
                    if not r["ID"]:
                        r["ID"] = "r_" + str(len(self.references) + 1)
                    self.references = pd.concat(
                            (self.references, pd.DataFrame.from_records([r])), 
                            ignore_index=True)
                    logging.debug("Added reference " + str(r))
            else:
                reference_list[i]["ID"] = ""
        return reference_list

    def get_reference(self, id):
        match = self.references[self.references["ID"] == id]
        if len(match) == 0:
            return None
        if len(match) > 1:
            logging.warn(f"Multiple matches for reference ID: {id}")
        return match.to_dict('records')[0]

    def find_reference(self, reference):
        for k in ["ID", "author", "title", "issued", "containerTitle", 
                "volume", "issue", "page", "citation"]:
            if k not in reference:
                reference[k] = ""

        match = self.references[
            (self.references["author"] == reference["author"]) &
            (self.references["title"] == reference["title"]) &
            (self.references["issued"] == reference["issued"]) &
            (self.references["containerTitle"] == reference["containerTitle"]) &
            (self.references["volume"] == reference["volume"]) &
            (self.references["issue"] == reference["issue"]) &
            (self.references["page"] == reference["page"]) &
            (self.references["citation"] == reference["citation"])]

        if len(match) == 0:
            return None
        return match.iloc[0]


    def start_name_bundle(self, accepted, incertae_sedis=False, sic=False):
        return NameBundle(self, accepted, incertae_sedis, sic)

    #
    # Ensure one or more names are included in names dataframe and update
    # taxa and synonyms dataframes as necessary
    # 
    # Parameters:
    #   name_dict - list of dictionaries - each dictionary contains 
    #       values keyed by terms from name_headings 
    #
    # Behaviour:
    #   Add names, taxon and synonyms to package and update bundle records
    #   with id values    
    #
    def add_names(self, bundle, parent):
        # Tag any names that are already known with the appropriate id
        self.prepare_bundle(bundle)

        # Ensure all names have IDs
        bundle.accepted = self.identify_name(bundle.accepted)
        if bundle.species:
            bundle.species = self.identify_name(bundle.species)
        if bundle.species_synonym:
            bundle.species_synonym = self.identify_name(bundle.species_synonym)
        for i in range(len(bundle.synonyms)):
            bundle.synonyms[i] = self.identify_name(bundle.synonyms[i])

        # Fix basionyms
        for i in range(len(bundle.synonyms)):
            if not bundle.synonyms[i]["basionymID"]:
                if not bundle.synonyms[i]["authorship"].startswith("("):
                    bundle.synonyms[i] = self.set_basionymid( 
                            bundle.synonyms[i],
                            bundle.synonyms[i]["ID"])
                    for j in range(len(bundle.synonyms)):
                        if i != j and not bundle.synonyms[j]["basionymID"] and \
                                self.same_basionym(
                                    bundle.synonyms[i], bundle.synonyms[j]):
                            bundle.synonyms[j] = self.set_basionymid( 
                                    bundle.synonyms[j],
                                    bundle.synonyms[i]["basionymID"])
                    if not bundle.synonyms[i]["basionymID"] and \
                            self.same_basionym(
                                bundle.synonyms[i], bundle.accepted):
                        # Will fix these below and map them to the basionymID
                        # when set.
                        bundle.synonyms[i] = self.set_basionymid( 
                                bundle.synonyms[i],
                                bundle.accepted["basionymID"])
                else:
                    # Hope we find the basionym later
                    pass

        if not bundle.accepted["basionymID"]:
            bundle.accepted = self.fix_basionymid(
                    bundle.accepted, bundle.synonyms)
        if bundle.species and not bundle.species["basionymID"]:
            bundle.species = self.fix_basionymid(bundle.species, bundle.synonyms)
        if bundle.species_synonym and not bundle.species_synonym["basionymID"]:
            bundle.species_synonym = self.set_basionymid(
                    bundle.species_synonym, bundle.species["basionymID"])

        if bundle.accepted["basionymID"] != bundle.accepted["ID"]:
            self.names.loc[self.names["basionymID"] == bundle.accepted["ID"], 
                    ["basionymID"]] = bundle.accepted["basionymID"]

        if bundle.species:
            taxon = self.add_taxon(bundle.species, parent)
            parent = taxon["ID"]
            bundle.species_taxon_id = parent
            if bundle.species_synonym:
                self.add_synonym(parent, bundle.species_synonym["ID"])

        taxon = self.add_taxon(bundle.accepted, parent, bundle.incertae_sedis)        
        bundle.accepted_taxon_id = taxon["ID"]
        for i in range(len(bundle.synonyms)):
            self.add_synonym(bundle.accepted_taxon_id, bundle.synonyms[i]["ID"])

        return

    def add_name_relation(self, name_relation):
        if self.name_relations is None:
            logging.info("Adding name_relation table")
            self.name_relations = pd.DataFrame(
                    columns=name_relation_headings)
        
        if "nameID" not in name_relation or \
                name_relation["nameID"] not in self.names["ID"].values:
            self.issue("Type material must be associated with a valid name ID")
            return None
        
        if "relatedNameID" not in name_relation or \
                name_relation["relatedNameID"] not in self.names["ID"].values:
            self.issue(
                "Type material must be associated with a valid related name ID")
            return None

        self.name_relations = \
                pd.concat((self.name_relations, 
                    pd.DataFrame.from_records([name_relation])), 
                            ignore_index=True)
        return name_relation

    def add_type_material(self, type_material):
        if self.type_materials is None:
            logging.info("Adding type_material table")
            self.type_materials = pd.DataFrame(
                    columns=type_material_headings)
        
        if "nameID" not in type_material or \
                type_material["nameID"] not in self.names["ID"].values:
            self.issue("Type material must be associated with a valid name ID")
            return None

        self.type_materials = \
                pd.concat((self.type_materials, 
                    pd.DataFrame.from_records([type_material])), ignore_index=True)
        return type_material

    def add_distribution(self, distribution):
        if self.distributions is None:
            logging.info("Adding distribution table")
            self.distributions = pd.DataFrame(
                    columns=distribution_headings)
        
        if "taxonID" not in distribution or \
                distribution["taxonID"] not in self.taxa["ID"].values:
            self.issue("Distribution must be associated with a valid taxon ID")
            return None

        self.distributions = \
                pd.concat((self.distributions, 
                    pd.DataFrame.from_records(
                            [distribution])), ignore_index=True)
        return distribution

    #
    # Add extra names
    def prepare_bundle(self, bundle):
        # Identify the species that should exist to allow the accepted
        # name to be grafted.
        # E.g. Aus bus cus --> Aus bus
        if self.species_from_trinomials:
            if bundle.accepted["infraspecificEpithet"]:
                bundle.species = \
                        bundle.derive_name(bundle.accepted, 
                            {"rank": "species",
                                "infraspecificEpithet": ""}, 
                            bundle.accepted_sic)
                # If this is a nominate subspecies, we can retain the 
                # authorship. Otherwise, check if we already know it.
                if bundle.accepted["infraspecificEpithet"] != \
                        bundle.accepted["specificEpithet"]:
                    bundle.species["authorship"] = ""
                    match = self.names[
                        (self.names["scientificName"] == \
                                bundle.species["scientificName"]) & \
                        (self.names["rank"] == "species")]
                    for i in range(len(match)):
                        if len(self.taxa[self.taxa["nameID"] == \
                                match.iloc[i]["ID"]]) > 0:
                            bundle.species["authorship"] = \
                                    match.iloc[i]["authorship"]
                            break

        # Add subspecies names for each infrasubspecific name.
        # E.g. Aus bus var. cus --> Aus bus cus
        if self.subspecies_from_trinomials:
            for i in range(len(bundle.synonyms)):
                if self.is_infrasubspecific(bundle.synonyms[i]):
                    bundle.add_synonym(
                            bundle.derive_name(bundle.synonyms[i], 
                                {"rank": "subspecies"}))

        # Add subgenus-free versions of names.
        # E.g. Aus (Aus) bus  --> Aus bus
        if self.subgenus_free_synonyms:
            if bundle.accepted["infragenericEpithet"]:
                bundle.add_synonym(
                        bundle.derive_name(bundle.accepted, 
                            {"infragenericEpithet": ""}, 
                            bundle.accepted_sic))

            for i in range(len(bundle.synonyms)):
                if bundle.synonyms[i]["infragenericEpithet"]:
                    bundle.add_synonym(
                            bundle.derive_name(bundle.synonyms[i], 
                                {"infragenericEpithet": ""}))

            if bundle.species and bundle.species["infragenericEpithet"]:
                bundle.species_synonym = \
                        bundle.derive_name(bundle.species, 
                            {"infragenericEpithet": ""})

    def add_synonym(self, taxon_id, name_id):
        match = self.synonyms[(self.synonyms["taxonID"] == taxon_id) & \
                (self.synonyms["nameID"] == name_id)]
        if len(match) > 0:
            return match.iloc[0]
        synonym_row = { 
            "taxonID": taxon_id, "nameID": name_id, "status": "synonym" }
        self.synonyms = pd.concat((self.synonyms, pd.DataFrame.from_records([synonym_row])), ignore_index = True)
        return

    def add_taxon(self, name, parent, incertae_sedis=False):
        match = self.taxa[self.taxa["nameID"] == name["ID"]]
        if len(match) > 0:
            if len(match) > 1:
                self.issue("More than one taxon matches name " + name["ID"])
            match = match.iloc[0]
            if match["uninomial"] != name["uninomial"]:
                self.issue("Uninomial mismatch (" + match["uninomial"] + \
                        " vs " + name["uninomial"] + ") for taxon id "
                        + match["ID"])
            if "species" in match and match["species"] and \
                    not name["scientificName"].startswith(match["species"]):
                self.issue("Species mismatch (" + match["species"] + \
                        " vs " + name["scientificName"] + ") for taxon id " 
                        + match["ID"])
            if match["parentID"] and parent and match["parentID"] != parent:
                match = self.taxa[self.taxa["ID"] == match["parentID"]].iloc[0]
                parentA = self.names[self.names["ID"] 
                            == match["nameID"]].iloc[0]
                match = self.taxa[self.taxa["ID"] == parent].iloc[0]
                parentB = self.names[self.names["ID"] 
                            == match["nameID"]].iloc[0]
                self.issue("Parent mismatch for " + name["rank"] + " " + \
                        name["scientificName"] + " " + name["authorship"] + \
                        ": was: " + parentA["scientificName"] + " " + \
                        parentA["authorship"] + ", now: " + \
                        parentB["scientificName"] + " " + parentB["authorship"])
            return match
        else:
            parent_rank = ""
            parent_name = ""
            if parent:
                match = self.taxa[self.taxa["ID"] == parent]
                if len(match) == 0:
                    self.issue("Parent taxon not found " + parent)
                    taxon_row = self.default_taxon
                else:
                    taxon_row = match.iloc[0].copy()
                    parent_name_row = \
                            self.names[self.names["ID"] == taxon_row["nameID"]]
                    parent_rank = parent_name_row.iloc[0]["rank"]
                    parent_name = parent_name_row.iloc[0]["scientificName"]
            else:
                taxon_row = self.default_taxon
            taxon_row["ID"] = str(len(self.taxa) + 1)
            taxon_row["parentID"] = parent if parent else ""
            taxon_row["nameID"] = name["ID"]
            if parent_rank in [ "kingdom", "phylum", "class", "order", 
                        "superfamily", "family", "subfamily" ]:
                taxon_row[parent_rank] = parent_name
            taxon_row["uninomial"] = name["uninomial"]
            if self.is_species_group(name):
                taxon_row["species"], _ = \
                        self.construct_species_rank_name(
                            name["genus"], name["infragenericEpithet"],
                            name["specificEpithet"], "", "")
                taxon_row["genus"] = name["genus"]
            else:
                taxon_row["species"] = ""
            if incertae_sedis:
                taxon_row["provisional"] = True
                if "remarks" in taxon_row and taxon_row["remarks"]:
                    taxon_row["remarks"] = \
                        f"Incertae sedis - {taxon_row['remarks']}"
                else:
                    taxon_row["remarks"] = "Incertae sedis"
            self.taxa = pd.concat((self.taxa, pd.DataFrame.from_records([taxon_row])), ignore_index=True)

            return taxon_row

    def modify_taxon(self, taxon_id, properties):
        self.taxa.loc[self.taxa["ID"] == taxon_id, list(properties.keys())] \
                = list(properties.values())

    def modify_name(self, name_id, properties):
        self.names.loc[self.names["ID"] == name_id, list(properties.keys())] \
                = list(properties.values())

    def identify_name(self, name):
        match = None
        if name["rank"] != "species" or not self.allow_repeated_binomials:
            match = self.find_name_record(name)
        if match is not None:
            name["ID"] = match["ID"]
            name["basionymID"] = match["basionymID"]
        else:
            if "ID" not in name or not name["ID"]:
                name["ID"] = str(len(self.names) + 1)
            self.names = pd.concat((self.names, pd.DataFrame.from_records([name])), ignore_index = True)
        return name

    def same_basionym(self, a, b):
        authorship_a = self.get_original_authorship(a["authorship"])
        authorship_b = self.get_original_authorship(b["authorship"])
        epithet_a = a["infraspecificEpithet"] if a["infraspecificEpithet"] \
                else a["specificEpithet"]
        epithet_b = b["infraspecificEpithet"] if b["infraspecificEpithet"] \
                else b["specificEpithet"]
        return authorship_a == authorship_b and \
                self.remove_gender(epithet_a) == self.remove_gender(epithet_b)

    def remove_gender(self, epithet):
        for suffix in [ "a", "us", "um", "is", "e", "os", "on" ]:
            if epithet.endswith(suffix):
                return epithet[0:len(epithet) - len(suffix)]
        return epithet

    def get_original_authorship(self, authorship):
        if authorship.startswith("(") and ")" in authorship:
            return authorship[1:authorship.index(")")]
        return authorship

    def epithet_and_authorship_match(self, name, epithet, authorship):
        epithet = self.remove_gender(epithet)
        return name["authorship"] == authorship and \
                (self.remove_gender(name["infraspecificEpithet"]) == epithet \
                    or (self.remove_gender(name["specificEpithet"]) == epithet \
                        and not name["infraspecificEpithet"]))

    def set_basionymid(self, name, basionymid):
        self.names.loc[self.names["ID"] == name["ID"], ["basionymID"]] = basionymid
        name["basionymID"] = basionymid
        return name

    def fix_basionymid(self, name, synonyms):
        original_authorship = self.get_original_authorship(name["authorship"])
        if name["authorship"] != original_authorship:
            epithet = name["infraspecificEpithet"] \
                    if name["infraspecificEpithet"] \
                    else name["specificEpithet"]
            for i in range(len(synonyms)):
                if self.epithet_and_authorship_match(
                        synonyms[i], epithet, original_authorship):
                    name = self.set_basionymid(name, 
                            synonyms[i]["basionymID"])
                    break
        else:
            name = self.set_basionymid(name, name["ID"])
        return name

    def find_name_record(self, name):
        return self.find_name( 
            name["scientificName"], name["authorship"], name["rank"])

    def get_name(self, id):
        match = self.names[self.names["ID"] == id]
        if len(match) == 0:
            return None
        if len(match) > 1:
            logging.warn(f"Multiple matches for name ID: {id}")
        return match.to_dict('records')[0]

    def find_name(self, scientificName, authorship, rank):
        if authorship is None:
            match = self.names[
                    (self.names["scientificName"] == scientificName) &
                    (self.names["rank"] == rank)]
        else:
            match = self.names[
                    (self.names["scientificName"] == scientificName) &
                    (self.names["authorship"] == authorship) &
                    (self.names["rank"] == rank)]
        if len(match) > 1:
            self.issue("Multiple name records for " + \
                            f"{rank} {scientificName} {str(authorship)}")
        if len(match) == 1:
            return match.iloc[0]
        return None

    def find_taxon(self, scientificName, authorship, rank):
        name = self.find_name(scientificName, authorship, rank)
        if name is not None:
            match = self.taxa[self.taxa["nameID"] == name["ID"]]
            if len(match) > 1:
                self.issue("Multiple taxa for name " + str(name))
            if len(match) > 0:
                return match.iloc[0]
        return None
    
    def find_names(self, properties, to_dict=False):
        for k in properties.keys():
            if k not in self.names.columns:
                self.issue(f"Unknown name property: {k}")
                return None
            
        names = self.names
        for k in properties.keys():
            names = names[names[k] == properties[k]]
            if len(names) == 0:
                return None
            
        if to_dict:
            return names.to_dict('records')
        return names

    def get_taxon(self, id):
        match = self.taxa[self.taxa["ID"] == id]
        if match.empty:
            return None
        if len(match.index  ) > 1:
            logging.warn(f"Multiple matches for taxon ID: {id}")
        return match.to_dict('records')[0]
    
    def get_synonyms(self, taxonID, to_dict=False):
        match = self.synonyms[self.synonyms["taxonID"] == taxonID]
        if match.empty:
            return None
        if to_dict:
            return match.to_dict('records')
        return match

    def get_synonymy(self, nameID, to_dict=False):
        match = self.taxa[self.taxa["nameID"] == nameID]
        if match.empty:
            match = self.synonyms[self.synonyms["nameID"] == nameID]
            if match.empty:
                return None, None
            if len(match.index) > 1:
                logging.warn(f"Multiple synonyms for name {nameID}")
            taxonID = match.iloc[0]["taxonID"]
            match = self.taxa[self.taxa["ID"] == taxonID]
            if match.empty:
                logging.error(f"Taxon not found for taxonID {taxonID}")
                return None, None
            if len(match.index) > 1:
                logging.warn(f"Multiple taxa found for taxonID {taxonID}")
            taxon = match.iloc[0]
        else:
            if len(match.index) > 1:
                logging.warn(f"Multiple taxa for name {nameID}")
            taxon = match.iloc[0]
        match = self.names[self.names["ID"] == taxon["nameID"]]
        if match.empty:
            logging.error(f"Name not found for nameID {taxon['nameID']}")
            accepted = None
        else:
            if len(match.index) > 1:
                logging.warn(f"Multiple matches for name ID: {id}")
            if to_dict:
                accepted = match.to_dict('records')[0]
            else:
                accepted = match.iloc[0]
        synonyms = self.get_synonyms(taxon['ID'])
        if synonyms.empty:
            synonymy = None
        else:
            synonymy = self.names[self.names["ID"].isin(synonyms["nameID"])]
            if to_dict:
                synonymy = synonymy.to_dict('records')

        return accepted, synonymy

    def construct_species_rank_name(self, g, sg, s, ss, marker):
        if not g:
            self.issue("Missing genus for species rank name " 
                                    + str([g, sg, s, ss, marker]))
            g = "?"
        if not s:
            self.issue("Missing specific epithet for species rank name " 
                                    + str([g, sg, s, ss, marker]))
            g = "?"

        if marker == "species":
            marker = None

        sn = [g, s]

        if ss:
            rank = "subspecies"
            sn.append(ss)
            if marker:
                if marker in ["subspecies", "subsp.", "sub.", "ss."]:
                    marker = "subsp." if self.code == "ICBN" else None
                    rank = "subspecies"
                elif marker in ["varietas", "variety", "var.", "v."]:
                    marker = "var."
                    rank = "variety"
                elif marker in ["subvarietas", "subvariety", "subvar."]:
                    marker = "subvar."
                    rank = "subvariety"
                elif marker in ["forma", "form", "f."]:
                    marker = "f."
                    rank = "form"
                elif marker in ["aberratio", "aberration", "ab."]:
                    marker = "ab."
                    rank = "aberration"
                elif marker in [ "infrasubspecific name" ]:
                    marker = ""
                    rank = "infrasubspecific name"
                else:
                    self.issue("Unknown rank marker " + marker)
                    rank = marker
                if marker:
                    sn.insert(2, marker)
                    marker = None
        else:
            rank = "species"

        if sg:
            if not sg.startswith("("):
                sg = "(" + sg + ")"
            sn.insert(1, sg)

        return " ".join(sn), rank

    def construct_authorship(self, a, y, is_basionym):
        if not a:
            self.issue("No author supplied")
            a = "?"
        if not y:
            self.issue("No year supplied")
            y = "?"

        year_parts = y.split()
        if len(year_parts) > 0:
            y = year_parts[0]
            year_publication = " ".join(year_parts[1:])
        else:
            year_publication = y
            
        if is_basionym:
            return a + ", " + y, y, year_publication
        else:
            return "(" + a + ", " + y + ")", y, year_publication

    def is_species_group(self, name):
        if name["specificEpithet"]:
            return True
        return False

    def is_infrasubspecific(self, name):
        return "rank" in name and \
                name["rank"] in ["variety", "form", "aberration"]

    #
    # Write dataframes as COLDP CSV files
    # 
    # Behaviour:
    #   Ensure that <folder>/<name>/ exists and write all dataframes as CSV.
    #
    def save(self, destination = None, name = None):
        if destination is None:
            destination = self.folder
        if destination is None or not os.path.exists(destination):
            logging.critical(
                f"Can't save package. Folder does not exist: {destination}")
            return
        logging.debug("Saving COLDP")    
        if name is None:
            name = self.name
        coldp_folder = os.path.join(destination, name)
        if not os.path.exists(coldp_folder):
            os.mkdir(coldp_folder)
        for item, value in { 
                "reference": self.references, 
                "name": self.names, 
                "taxon": self.taxa,
                "synonym": self.synonyms,
                "typematerial": self.type_materials, 
                "distribution": self.distributions, 
                "namerelation": self.name_relations,
                "speciesinteraction": self.species_interactions,
                "issue": self.issues 
                }.items():
            if value is not None and len(value) > 0:
                value.replace('', np.nan, inplace=True)
                value.dropna(how='all', axis=1, inplace=True)
                value.replace(np.nan, '', inplace=True)
                file_name = os.path.join(coldp_folder, item + ".csv")
                value.to_csv(file_name, index=False)
        logging.info("COLDP saved to " + coldp_folder)
        return

    def issue(self, message, **record):
        if self.issues is None:
            self.issues = pd.DataFrame(columns=issues_headings)

        context = self.context if self.context else ""
        record = {"issue": message, "context": context }
        
        if self.issues_to_stdout:
            logging.error(context + ": " + message)
        
        self.issues = pd.concat((self.issues, pd.DataFrame.from_records([record])), ignore_index=True)
            

        
if __name__ == "__main__":

    logging.basicConfig(filename='COLDP.log', level=logging.DEBUG)

    default_taxon = {
        "kingdom": "Animalia",
        "phylum": "Arthropoda",
        "class": "Insecta",
        "order": "Lepidoptera",
        "scrutinizer": "Geometridae Working Group",
        "provisional": "false",
        "extinct": "false",
        "lifezone": "terrestrial",
        "temporalRangeEnd": "Holocene"
    }

    for f in os.listdir("./test/output/"):
        shutil.rmtree(os.path.join("./test/output/", f))

    coldp = COLDP("./test/output/", "fake", 
                    default_taxon_record=default_taxon,
                    insert_species_for_trinomials=True,
                    create_subspecies_for_infrasubspecifics=True,
                    create_synonyms_without_subgenus=True,
                    issues_to_stdout=True)

    coldp.add_references([
        {"author": "White, T.H.", "issued": "1950", "title": "The Once and Future King", "containerTitle": "Penguin Classics", "volume": "34", "page": "1-756"},
        {"author": "Gielis, C.", "issued": "2022", "title": "Moths of Bhutan", "page": "1-382"}
    ])

    coldp.add_references([
        {"author": "Gielis, C.", "issued": "2022", "title": "Moths of Bhutan", "page": "1-382"},
        {"author": "Gielis, C.", "issued": "2022", "title": "More Moths of Bhutan", "page": "1-112"}
    ])

    nb = coldp.start_name_bundle({
        "genus": "Aus",
        "infragenericEpithet": "Xus",
        "specificEpithet": "bus",
        "infraspecificEpithet": "bus",
        "authorship": "(Jones, 1952)",
        "code": "ICZN",
        "status": "established"
    })

    nb.add_synonym({
        "genus": "Xus",
        "specificEpithet": "cus",
        "authorship": "Jones, 1952",
        "code": "ICZN",
        "status": "established"
    })

    nb.add_synonym({
        "genus": "Aus",
        "specificEpithet": "cus",
        "authorship": "(Jones, 1952)",
        "code": "ICZN",
        "status": "established"
    })

    nb.add_synonym({
        "genus": "Yus",
        "specificEpithet": "bus",
        "authorship": "Jones, 1952",
        "code": "ICZN",
        "status": "established"
    })

    coldp.add_names(nb, None)

    coldp.modify_taxon(nb.accepted_taxon_id, {"extinct": "true", "lifezone": "marine"})

    coldp.save()

    coldp = COLDP("./test/input", "Tonzidae", 
                basionyms_from_synonyms=True,
                classification_from_parents=True,
                issues_to_stdout=True)

    coldp.save("./test/output")

    coldp = COLDP("../../Dropbox/NetBeans Projects/Tineidae", "NHM_Tineidae_2011-11-21", 
                basionyms_from_synonyms=True,
                classification_from_parents=True,
                issues_to_stdout=True)

    coldp.save("./test/output")
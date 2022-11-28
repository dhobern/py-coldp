#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#----------------------------------------------------------------------------
# Created By  : Donald Hobern, dhobern@gmail.com
# Created Date: 2022-11-27
# version ='1.0'
# ---------------------------------------------------------------------------
""" 
COLDP.py

COLDP is a Python class for creating and serialising Catalogue of Life
Data Package (COLDP) files.

COLDP is the preferred format within Catalogue of Life for organising taxonomic 
checklist datasets. See: https://github.com/CatalogueOfLife/coldp 
"""
# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import pandas
import logging
import os
import re

reference_headings = [ 
    "ID", "author", "title", "issued", "containerTitle", "volume", "issue", 
    "page", "link", "citation" ]

name_headings = [ 
    "ID", "basionymID", "scientificName", "authorship", "rank", "uninomial", 
    "genus", "infragenericEpithet", "specificEpithet", "infraspecificEpithet", 
    "referenceID", "publishedInPage", "publishedInYear", "code", "status", 
    "remarks", "link" ]

taxon_headings = [ 
    "ID", "parentID", "nameID", "scrutinizer", "scrutinizerDate", "provisional",
     "referenceID", "extinct", "temporalRangeEnd", "lifezone", "kingdom", 
     "phylum", "class", "order", "superfamily", "family", "subfamily", "tribe", 
     "genus", "uninomial", "species", "remarks" ]

synonym_headings = [ 
    "taxonID", "nameID", "accordingToID", "status", "referenceID", "remarks" ]

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

issues_headings = [ "context", "issue" ]

uninomial_pattern = re.compile('^[A-Z][a-z]+$')
epithet_pattern = re.compile('^[a-z]-?[a-z]+$')
name_pattern = re.compile('^[A-Za-z]-?[a-z]+$')

#
# CLASS: NameBundle
#
# Wrapper to manage set of associated names, normally an original and an
# accepted name. The bundle handles the logic of associate name variations.
# 
# Usage:
#   The minimal usage is to create a new bundle with an accepted name. One 
#   or more synonyms can also be supplied using the add_synonym() method.
#   The prepare() method augments the bundle with other associated names to
#   support intelligent handling by the COLDP package. 
#   (variant) 
#
class NameBundle:

    def __init__(self, coldp, accepted, incertae_sedis=False):
        self.coldp = coldp
        if "rank" not in accepted:
            self.coldp.issue(
                "Accepted name missing rank, assume species: " + str(accepted))
            accepted["rank"] = "species"
        self.accepted = self.normalise_name(accepted)
        self.incertae_sedis = incertae_sedis
        self.accepted_taxon_id = None
        self.synonyms = []
        self.species = None
        self.species_synonym = None
        self.accepted_species_id = None

    def add_synonym(self, synonym, sic=False):
        self.synonyms.append(self.normalise_name(synonym, sic))

    def normalise_name(self, name, sic=True):
        for k in [
                "uninomial", "genus", "infragenericEpithet", "specificEpithet", 
                "infraspecificEpithet", "authorship", "rank", "basionymID" ]:
            if k not in name:
                name[k] = ""
        if not sic:
            for k in [ "uninomial", "genus", "infragenericEpithet" ]:
                if len(name[k]) > 0:
                    if not name_pattern.match(name[k]):
                        self.coldp.issue(
                            "Invalid pattern for supraspecific name: \"" + name[k] + "\"")
                    elif name[k][0].islower():
                        self.coldp.issue(
                            "Lowercase initial letter for supraspecific name: " 
                                + name[k])
                        name[k] = name[k][0].upper() + name[k][1:]
            for k in [ "specificEpithet", "infraspecificEpithet" ]:
                if len(name[k]) > 0:
                    if not name_pattern.match(name[k]):
                        self.coldp.issue("Invalid pattern for epithet: " + name[k])
                    elif name[k][0].isupper():
                        self.coldp.issue(
                            "Uppercase initial letter for epithet: " + name[k])
                        name[k] = name[k][0].lower() + name[k][1:]
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
        return name

    def derive_name(self, name, values):
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
        return self.normalise_name(values)

#
# CLASS: COLDP
#
#   Wrapper for set of Pandas dataframes for CSV tables in Catalogue of Life
#   Data Package.
#
#   Constructor:
#       folder - name of folder in which to save COLDP package
#       name - name for COLDP package (files will be saved as 
#           <folder>/<name>/*.csv)
#
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
        self.issues_to_stdout = False
        if kwargs:
            self.set_options(kwargs)
        self.folder = folder
        self.name = name
        self.names = pandas.DataFrame(columns=name_headings)
        self.taxa = pandas.DataFrame(columns=taxon_headings)
        self.synonyms = pandas.DataFrame(columns=synonym_headings)
        self.references = pandas.DataFrame(columns=reference_headings)
        self.type_materials = None
        self.distributions = None
        self.name_relations = None
        self.species_interactions = None
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
            elif key == "code" and value in ["ICZN", "ICBN"]:
                self.code = value

    def set_context(self, context):
        self.context = context

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
            r = reference_list[i]
            for k in ["author", "title", "issued", "containerTitle", "volume", 
                    "issue", "page", "citation"]:
                if k not in r:
                    r[k] = ""

            match = self.references[
                (self.references["author"] == r["author"]) &
                (self.references["title"] == r["title"]) &
                (self.references["issued"] == r["issued"]) &
                (self.references["containerTitle"] == r["containerTitle"]) &
                (self.references["volume"] == r["volume"]) &
                (self.references["issue"] == r["issue"]) &
                (self.references["page"] == r["page"]) &
                (self.references["citation"] == r["citation"])]
            if len(match) > 1:
                logging.critical(
                    "Multiple reference records exist for " + str(r))
            elif len(match) == 1:
                reference_list[i] = match.iloc[0]
                logging.debug(
                    "Matched reference ID " + reference_list[i]["ID"])
            else:
                if not r["ID"]:
                    r["ID"] = str(len(self.references) + 1)
                self.references = pandas.concat((self.references, pandas.DataFrame.from_records([r])), ignore_index=True)
                logging.debug("Added reference " + str(r))
        return reference_list

    def start_name_bundle(self, accepted, incertae_sedis=False):
        return NameBundle(self, accepted)

    #
    # Ensure one or more names are included in names dataframe and update
    # taxa and synonyms dataframes as necessary
    # 
    # Parameters:
    #   name_dict - list of dictionaries - each dictionary contains 
    #       values keyed by terms from name_headings 
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
                    for j in range(i+1, len(bundle.synonyms)):
                        if not bundle.synonyms[j]["basionymID"] and \
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

    def add_type_material(self, type_material):
        if self.type_materials is None:
            logging.info("Adding type_material table")
            self.type_materials = pandas.DataFrame(
                    columns=type_material_headings)
        
        if "nameID" not in type_material or \
                type_material["nameID"] not in self.names["ID"].values:
            self.issue("Type material must be associated with a valid name ID")
            return None

        self.type_materials = \
                pandas.concat((self.type_materials, 
                    pandas.DataFrame.from_records([type_material])), ignore_index=True)
        return type_material

    def add_distribution(self, distribution):
        if self.distributions is None:
            logging.info("Adding distribution table")
            self.distributions = pandas.DataFrame(
                    columns=distribution_headings)
        
        if "taxonID" not in distribution or \
                distribution["taxonID"] not in self.taxa["ID"].values:
            self.issue("Distribution must be associated with a valid taxon ID")
            return None

        self.distributions = \
                pandas.concat((self.distributions, 
                    pandas.DataFrame.from_records(
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
                                "infraspecificEpithet": ""})
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
                            {"infragenericEpithet": ""}))

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
        self.synonyms = pandas.concat((self.synonyms, pandas.DataFrame.from_records([synonym_row])), ignore_index = True)
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
            if match["species"] and \
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
                if taxon_row["remarks"]:
                    taxon_row["remarks"] = "Incertae sedis - " + taxon_row["remarks"]
                else:
                    taxon_row["remarks"] = taxon_row["remarks"]
                taxon_row["remarks"] = "Incertae sedis"
            self.taxa = pandas.concat((self.taxa, pandas.DataFrame.from_records([taxon_row])), ignore_index=True)

            return taxon_row

    def modify_taxon(self, taxon_id, properties):
        self.taxa.loc[self.taxa["ID"] == taxon_id, list(properties.keys())] \
                = list(properties.values())

    def identify_name(self, name):
        match = self.find_name_record(name)
        if match is not None:
            name["ID"] = match["ID"]
            name["basionymID"] = match["basionymID"]
        else:
            if "ID" not in name or not name["ID"]:
                name["ID"] = str(len(self.names) + 1)
            self.names = pandas.concat((self.names, pandas.DataFrame.from_records([name])), ignore_index = True)
        return name

    def same_basionym(self, a, b):
        authorship_a = a["authorship"]
        if authorship_a.startswith("("):
            authorship_a = authorship_a[1: len(authorship_a) - 1]
        authorship_b = b["authorship"]
        if authorship_b.startswith("("):
            authorship_b = authorship_b[1: len(authorship_b) - 1]
        epithet_a = a["infraspecificEpithet"] if a["infraspecificEpithet"] \
                else a["specificEpithet"]
        epithet_b = b["infraspecificEpithet"] if b["infraspecificEpithet"] \
                else b["specificEpithet"]
        return authorship_a == authorship_b and epithet_a == epithet_b

    def epithet_and_authorship_match(self, name, epithet, authorship):
        return name["authorship"] == authorship and \
                (name["infraspecificEpithet"] == epithet or \
                    (name["specificEpithet"] == epithet and \
                        not name["infraspecificEpithet"]))

    def set_basionymid(self, name, basionymid):
        self.names.loc[self.names["ID"] == name["ID"], ["basionymID"]] = basionymid
        name["basionymID"] = basionymid
        return name

    def fix_basionymid(self, name, synonyms):
        if name["authorship"].startswith("("):
            bare_authorship = name["authorship"] \
                    [1:len(name["authorship"]) - 1]
            epithet = name["infraspecificEpithet"] \
                    if name["infraspecificEpithet"] else name["specificEpithet"]
            for i in range(len(synonyms)):
                if self.epithet_and_authorship_match(
                        synonyms[i], epithet, bare_authorship):
                    name = self.set_basionymid(
                            name, synonyms[i]["basionymID"])
        else:
            name = self.set_basionymid(name, name["ID"])
        return name

    def find_name_record(self, name):
        return self.find_name( 
            name["scientificName"], name["authorship"], name["rank"])

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
            self.issue("Multiple name records for " + rank + " " + \
                    scientificName + " " + authorship)
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
    def save(self):
        if not os.path.exists(self.folder):
            logging.critical(
                "Can't save package. Folder does not exist: " + self.folder)
            return
        logging.debug("Saving COLDP")    
        coldp_folder = os.path.join(self.folder, self.name)
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
            if value is not None:
                file_name = os.path.join(coldp_folder, item + ".csv")
                value.to_csv(file_name, index=False)
        logging.info("COLDP saved to " + coldp_folder)
        return

    def issue(self, message, **record):
        if self.issues is None:
            self.issues = pandas.DataFrame(columns=issues_headings)

        context = self.context if self.context else ""
        record = {"issue": message, "context": context }
        
        if self.issues_to_stdout:
            logging.error(context + ": " + message)
        
        self.issues = pandas.concat((self.issues, pandas.DataFrame.from_records([record])), ignore_index=True)
            

        
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

    coldp = COLDP("./", "coldp", 
                    default_taxon_record=default_taxon,
                    insert_species_for_trinomials=True,
                    create_subspecies_for_infrasubspecifics=True,
                    create_synonyms_without_subgenus=True)

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
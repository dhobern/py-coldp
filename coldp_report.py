from coldp import COLDP
import csv

ranks = {
    "kingdom": "01",
    "phylum": "02",
    "class": "03",
    "order": "04",
    "suborder": "05",
    "superfamily": "06",
    "family": "07",
    "subfamily": "08",
    "tribe": "09",
    "subtribe": "10",
    "genus": "11",
    "subgenus": "12",
    "species": "13",
    "subspecies": "14"
}

headings = [
    "family", "code", "rank", "taxonID", "taxonModified", "nameID", "scientificName", "authorship", "nameModified", "referenceID", "citation", "referenceModified"
]

def rank_position(rank):
    if rank in ranks:
        return ranks[rank]
    return "99"

def report_taxon(writer, coldp, taxon, family):
    name = coldp.get_name(taxon["nameID"])
    if name["rank"] == "family":
        family = name["scientificName"]
    if name["referenceID"]:
        reference = coldp.get_reference(name["referenceID"])
    else:
        reference = None
    if reference is not None:
        referenceID = reference["ID"]
        citation = reference["citation"]
    else:
        referenceID = ""
        citation = ""

    row = [family, "A", name["rank"], taxon["ID"], taxon["modified"], name["ID"], name["scientificName"], name["authorship"], name["modified"], referenceID, citation]
    writer.writerow(row)

    synonyms = coldp.get_synonyms(taxon["ID"], True)
    if synonyms is not None:
        names = []
        for synonym in synonyms:
            names.append(coldp.get_name(synonym["nameID"]))
        names = sorted(names, key=lambda x: (x["publishedInYear"] if x["publishedInYear"] else "9999") + x["scientificName"])
        for other_name in names:
            code = "B" if name["ID"] == name["basionymID"] else "C"
        if name["referenceID"]:
            reference = coldp.get_reference(other_name["referenceID"])
        else:
            reference = None
        if reference is not None:
            referenceID = reference["ID"]
            citation = reference["citation"]
        else:
            referenceID = ""
            citation = ""
        row = [family, code, name["rank"], taxon["ID"], taxon["modified"], name["ID"], other_name["scientificName"], other_name["authorship"], other_name["modified"], referenceID, citation]
        writer.writerow(row)
    children = coldp.get_children(taxon["ID"], True)
    if children is not None:
        print(f"{name['scientificName']} ({name['rank']})")
        for child in children:
            name = coldp.get_name(child["nameID"])
            child["position"] = rank_position(name["rank"]) + name["scientificName"]
        children = sorted(children, key=lambda x: x["position"])
        for child in sorted(children, key=lambda x: x["position"]):
            report_taxon(writer, coldp, child, family)

with open("report.csv", 'w', newline='', encoding='utf8') as f:
    writer= csv.writer(f, delimiter=',')
    writer.writerow(headings)
    coldp = COLDP("D:/stang/Downloads", "latest", issues_to_stdout=True)
    taxon = coldp.find_taxon("Lepidoptera", None, "order")
    report_taxon(writer, coldp, taxon, "")
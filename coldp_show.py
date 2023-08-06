from coldp import COLDP

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


def rank_position(rank):
    if rank in ranks:
        return ranks[rank]
    return "99"

def show_taxon(f, coldp, taxon, indent):
    name = coldp.get_name(taxon["nameID"])
    f.write(f"{indent}{name['scientificName']} {name['authorship']} ({name['rank']})\n")
    synonyms = coldp.get_synonyms(taxon["ID"], True)
    if synonyms is not None:
        names = []
        for synonym in synonyms:
            names.append(coldp.get_name(synonym["nameID"]))
        names = sorted(names, key=lambda x: (x["publishedInYear"] if x["publishedInYear"] else "9999") + x["scientificName"])
        for other_name in names:
            f.write(f"{indent}  *{other_name['scientificName']} {other_name['authorship']} ({other_name['rank']})\n")
    children = coldp.get_children(taxon["ID"], True)
    if children is not None:
        print(f"{indent}{name['scientificName']} ({name['rank']})")
        for child in children:
            name = coldp.get_name(child["nameID"])
            child["position"] = rank_position(name["rank"]) + name["scientificName"]
        children = sorted(children, key=lambda x: x["position"])
        for child in sorted(children, key=lambda x: x["position"]):
            show_taxon(f, coldp, child, indent + "  ")

with open("show.txt", 'w', encoding='utf8') as f:
    coldp = COLDP("D:/stang/Downloads", "latest", issues_to_stdout=True)
    taxon = coldp.find_taxon("Lepidoptera", None, "order")
    indent = ""
    show_taxon(f, coldp, taxon, indent)
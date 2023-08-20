from coldp import COLDP
import sys

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
    "subspecies": "14",
}


def rank_position(rank):
    if rank in ranks:
        return ranks[rank]
    return "99"


def get_extras(name, extras, separator=", "):
    if extras is None or len(extras) == 0:
        return ""
    return f" ({separator.join([name[x] for x in extras])})"


def show_taxon(f, coldp, taxon, indent, extras=["rank"], separator=", "):
    name = coldp.get_name(taxon["nameID"])
    f.write(
        f"{indent}{name['scientificName']} {name['authorship']}{get_extras(name, extras, separator)}\n"
    )
    synonyms = coldp.get_synonyms(taxon["ID"], True)
    if synonyms is not None:
        names = []
        for synonym in synonyms:
            names.append(coldp.get_name(synonym["nameID"]))
        names = sorted(
            names,
            key=lambda x: (x["publishedInYear"] if x["publishedInYear"] else "9999")
            + x["scientificName"],
        )
        for other_name in names:
            f.write(
                f"{indent}  *{other_name['scientificName']} {other_name['authorship']}{get_extras(other_name, extras, separator)}\n"
            )
    children = coldp.get_children(taxon["ID"], True)
    if children is not None:
        print(f"{indent}{name['scientificName']}{get_extras(name, extras, separator)}")
        for child in children:
            name = coldp.get_name(child["nameID"])
            child["position"] = rank_position(name["rank"]) + name["scientificName"]
        children = sorted(children, key=lambda x: x["position"])
        for child in sorted(children, key=lambda x: x["position"]):
            show_taxon(f, coldp, child, indent + "  ", extras, separator)


if len(sys.argv) > 2:
    with open("show.txt", "w", encoding="utf8") as f:
        coldp = COLDP(sys.argv[1], sys.argv[2], issues_to_stdout=True)
        taxon = coldp.find_taxon("Lepidoptera", None, "order")
        indent = ""
        extras = []
        for e in ["rank", "modified", "modifiedBy", "referenceID"]:
            if e in coldp.names.columns:
                extras.append(e)
        show_taxon(f, coldp, taxon, indent, extras, separator=",")
else:
    print("Usage: python coldp_show.py folder coldp_name")

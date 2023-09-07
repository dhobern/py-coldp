from coldp import COLDP
import sys
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
    "subspecies": "14",
}


def get_rank_position(rank):
    if rank in ranks:
        return ranks[rank]
    return "99"


def report_taxon(coldp, taxon, extras=[], log_level=None, list=None):
    if list is None:
        list = []

    name = coldp.get_name(taxon["nameID"])
    rank_position = get_rank_position(name["rank"])
    if log_level is not None and log_level >= rank_position:
        print(
            f"{' ' * (int(rank_position) - 1)}{name['scientificName']} ({name['rank']})"
        )
    row = [
        name["rank"],
        rank_position,
        "accepted",
        name["scientificName"],
        name["authorship"],
    ]
    row += [name[x] for x in extras]
    list.append(row)

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
            row = [
                other_name["rank"],
                get_rank_position(other_name["rank"]),
                "synonym",
                other_name["scientificName"],
                other_name["authorship"],
            ]
            row += [other_name[x] for x in extras]
            list.append(row)

    children = coldp.get_children(taxon["ID"], True)
    if children is not None:
        for child in children:
            name = coldp.get_name(child["nameID"])
            child["position"] = get_rank_position(name["rank"]) + name["scientificName"]
        children = sorted(children, key=lambda x: x["position"])
        for child in sorted(children, key=lambda x: x["position"]):
            report_taxon(coldp, child, extras=extras, log_level=log_level, list=list)

    return list


if __name__ == "__main__":
    if len(sys.argv) > 2:
        if len(sys.argv) > 4:
            name = sys.argv[3]
            rank = sys.argv[4]
        else:
            name = "Lepidoptera"
            rank = "order"
        with open("show.csv", "w", encoding="utf8", newline="") as f:
            writer = csv.writer(f, delimiter=",")

            coldp = COLDP(sys.argv[1], sys.argv[2], issues_to_stdout=True)
            taxon = coldp.find_taxon(name, None, rank)
            extras = []
            for e in ["modified", "modifiedBy", "referenceID"]:
                if e in coldp.names.columns:
                    extras.append(e)
            hierarchy = report_taxon(coldp, taxon, extras=extras, log_level="07")
            headings = [
                "rank",
                "rankLevel",
                "status",
                "scientificName",
                "authorship",
            ] + extras
            writer.writerow(headings)
            for row in hierarchy:
                writer.writerow(row)
    else:
        print(f"Usage: python {sys.argv[0]} folder coldp_name taxon rank")

from coldp import COLDP

c = COLDP("Test", insert_species_for_trinomials=True)

nb = c.start_name_bundle(
    {"rank": "genus", "uninomial": "Aus", "authorship": "Linnaeus, 1758"}
)
nb.add_synonym({"rank": "genus", "uninomial": "Bus", "authorship": "Meyrick, 1912"})
c.add_names(nb)

gid = nb.accepted_taxon_id

nb = c.start_name_bundle(
    {
        "rank": "subspecies",
        "genus": "Aus",
        "specificEpithet": "bus",
        "infraspecificEpithet": "bus",
        "authorship": "Linnaeus, 1758",
    }
)
c.add_names(nb, gid)

nb = c.start_name_bundle(
    {
        "rank": "subspecies",
        "genus": "Aus",
        "specificEpithet": "bus",
        "infraspecificEpithet": "cus",
        "authorship": "Linnaeus, 1758",
    }
)
c.add_names(nb, gid)

print(c.get_text_tree(gid, "      "))

c.save()

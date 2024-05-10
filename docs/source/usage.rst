Usage
=====

The following example illustrates basic usage. It creates a new COLDP
instance for the monotypic Lepidoptera family Tridentaformidae, including
name and taxon records for the family, genus and species, the original 
combination for the species, three distribution records for the species, and
references for all elements.

Running the same code a second time leaves the data unchanged since it 
locates and validates the existing records.

In this example, each reference is associated with an ID string provided
when the records are created, but all name and taxon ID strings are 
generated and managed by the COLDP instance.

Adding the synonym to the NameBundle for the species enables the COLDP instance
to create the relevant basionymID reference to the current combination.

The example names include a mixture with parsed name elements or a single
scientificName element. The final names table includes both formats for all
names.

The addition of the distribution records shows the pattern for adding other
COLDP classes that reference name (namerelation, typematerial) and taxon 
(distribution, speciesinteraction) records.

.. code-block:: python
        
    # Import COLDP class
    from coldp import COLDP

    # Default properties for all COLDP Taxon records created by COLDP instance
    taxon_defaults = {
        "kingdom": "Animalia",
        "phylum": "Arthropoda",
        "class": "Insecta",
        "order": "Lepidoptera",
        "status": "established",
    }

    # Create new COLDP instance with name Tridentaformidae
    #
    # This would load an existing COLDP instance from a folder named
    # Tridentaformidae in the current folder if it already exists
    coldp = COLDP("Tridentaformidae", default_taxon_record = taxon_defaults)

    # Add four COLDP Reference objects
    references = coldp.add_references([
        {
            "ID": "Braun_1923",
            "author": "Braun, A.F.",
            "issued": "1923",
            "title": "Microlepidoptera: Notes and New Species",
            "containerTitle": "Transactions of the American Entomological Society",
            "volume": "49",
            "issue": "2", 
            "page": "115-127",
            "link": "https://www.jstor.org/stable/25077087",
        },
        {
            "ID": "Davis_1978",
            "author": "Davis, D.R.",
            "issued": "1978",
            "title": "Two new genera of North American incurvariine moths (Lepidoptera: Incurvariidae)",
            "containerTitle": "The Pan-Pacific entomologist",
            "volume": "54",
            "issue": "2", 
            "page": "147-153",
            "link": "https://www.biodiversitylibrary.org/page/56100973",
        },
        {
            "ID": "Pohl_et_al_2019",
            "author": "Pohl, G.R., Landry, J.-F., Schmidt, B.C. & deWaard, J.R.",
            "issued": "2019",
            "title": "Lepidoptera of Canada",
            "containerTitle": "ZooKeys",
            "volume": "819",
            "page": "463-505",
            "link": "https://doi.org/10.3897/zookeys.819.27259",
        },
        {
            "ID": "Regier_et_al_2014",
            "author": "Regier, J.C., Mitter, C., Davis, D.R., Harrison, T.L., Sohn, J.-C., Cummings, M.P., Zwick, A. & Mitter, K.T.",
            "issued": "2015",
            "title": "A molecular phylogeny for the oldest (nonditrysian) lineages of extant Lepidoptera, with implications for classification, comparative morphology and life-history evolution",
            "containerTitle": "Systematic Entomology",
            "volume": "40",
            "issue": "4", 
            "page": "671–704",
            "link": "https://doi.org/10.1111/syen.12129",
        },
    ])

    # Add COLDP Name and Taxon records for family and get the ID string for the 
    # family Taxon record
    #
    # Name provided as a uninomial
    bundle = coldp.start_name_bundle({
        "rank": "family",
        "uninomial": "Tridentaformidae",
        "authorship": "Davis, 2014",
        "referenceID": "Regier_et_al_2014",
        "publishedInPage": "697",
    })
    coldp.add_names(bundle)
    family_id = bundle.accepted_taxon_id

    # Add COLDP Name and Taxon records for genus as child of the family and get 
    # the ID string for the genus Taxon record
    #
    # Name provided as a scientificName
    bundle = coldp.start_name_bundle({
        "rank": "genus",
        "scientificName": "Tridentaforma",
        "authorship": "Davis, 1978",
        "referenceID": "Davis_1978",
        "publishedInPage": "150",
    })
    coldp.add_names(bundle, family_id)
    genus_id = bundle.accepted_taxon_id

    # Add COLDP Name and Taxon records for species as child of the genus
    # with original combination as a synonym and get the ID string for the 
    # species Taxon record
    #
    # Accepted name provided as parsed elements. Synonym provided only as 
    # scientificName.
    bundle = coldp.start_name_bundle({
        "rank": "species",
        "genus": "Tridentaforma",
        "specificEpithet": "fuscoleuca",
        "authorship": "(Braun, 1923)",
        "referenceID": "Davis_1978",
        "publishedInPage": "150",
        "publishedInYear": "1978",
    })
    bundle.add_synonym({
        "rank": "species",
        "scientificName": "Lampronia fuscoleuca",
        "authorship": "Braun, 1923",
        "referenceID": "Braun_1923",
        "publishedInPage": "127",
    })
    coldp.add_names(bundle, genus_id)
    species_id = bundle.accepted_taxon_id

    # Add three distribution records for the species, each with a reference
    for area, referenceID in {"US-CA": "Braun_1923", "CA-AB": "Pohl_et_al_2019", "CA-BC": "Pohl_et_al_2019"}.items():
        distribution = coldp.add_distribution({
            "taxonID": species_id,
            "area": area,
            "gazetteer": "iso",
            "status": "native",
            "referenceID": referenceID,
        })

    # Save the COLDP instance to a Tridentaformidae subfolder in the current 
    # folder
    coldp.save()

    # Load the COLDP instance from the curent folder
    coldp = COLDP("Tridentaformidae")

    # Display the classification as a text tee
    print(coldp.get_text_tree(family_id))

    # Show content of dataframes
    print(coldp.references)
    print(coldp.names)
    print(coldp.taxa)
    print(coldp.synonyms)
    print(coldp.distributions)

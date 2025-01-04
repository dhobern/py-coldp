# py-coldp
 Python tools for working with taxonomic checklists organised as Catalogue of Life Data Package (COLDP) format

# Overview
py-coldp is a Python package to facilitate creation, manipulation, editing and serialisation of taxonomic checklists in the [Catalogue of Life Data Package](https://github.com/CatalogueOfLife/coldp ) format.

The package includes two classes:
* **COLDP** - A COLDP package loaded as a set of Pandas dataframes
* **NameBundle** - A helper class to simplify addition of taxon names with sets of associated synonyms to a COLDP instance

# Class: COLDP
The main **COLDP** class instantiates a COLDP package in memory as a set of Pandas dataframes. An instance may be initialised from the contents of a folder containing a set of COLDP-compliant CSV or tab-delimited data files or alternatively can be initialised as an empty instance in memory. The class includes many methods for inserting new data, editing existing records and querying the contents of the package. The instance can then be saved as a set of CSV files in a named folder .

# Class: NameBundle

The **NameBundle** class brings together a scientific name and its synonyms so that these can be added together and the COLDP package can automatically manage their relationships NameBundle objects are normally created using the COLDP.**start_name_bundle()** method. 

At minimum a NameBundle is initialised with a dictionary holding a set of COLDP name record values. The scientific name represented by this dictionary should be the accepted name for a species or other taxon. Synonyms may then be added to the NameBundle. 

Once all names are included, the NameBundle can be added to the COLDP object via the COLDP.**add_names()** method. This adds name, taxon and synonym records for the set of names supplied. COLDP options may expand the set of added synonyms to include variant formats or may trigger the addition of one or taxa that are implicit in the accepted name. 

The COLDP object will automatically manage record identifiers and the basionymID for any name records that are combinations of another name in the set.

# Installation
```console
pip install py-coldp
```

# Documentation
[Browsable documentation for COLDP version: 2024.12.3](https://html-preview.github.io/?url=https://github.com/dhobern/py-coldp/blob/main/docs/build/html/index.html)

[PDF documentation for COLDP version: 2024.12.3](https://github.com/dhobern/py-coldp/blob/main/docs/build/simplepdf/COLDP.pdf)

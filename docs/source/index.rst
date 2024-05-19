.. Root for documentation for the coldp package

coldp
=====
Python tools for working with taxonomic checklists organised as Catalogue of Life Data Package (COLDP) format

Overview
--------
coldp is a Python package to facilitate creation, manipulation, editing and serialisation of taxonomic checklists in the `COL Data Package (COLDP) <https://github.com/CatalogueOfLife/coldp>`_ format.

The package includes three classes:

* **coldp.COLDP** - A COLDP package loaded as a set of Pandas dataframes
* **coldp.NameBundle** - A helper class to simplify addition of taxon names with sets of associated synonyms to a COLDP instance
* **coldp.IdentifierPolicy** - An internal class to manage ID values in COLDP dataframes

The main **COLDP** class instantiates a COLDP package in memory as a set of Pandas dataframes. An instance may be initialised from the contents of a folder containing a set of COLDP-compliant CSV or tab-delimited data files or alternatively can be initialised as an empty instance in memory. The class includes many methods for inserting new data, editing existing records and querying the contents of the package. The instance can then be saved as a set of CSV files in a named folder .

The **NameBundle** class brings together a scientific name and its synonyms so that these can be added together and the COLDP package can automatically manage their relationships NameBundle objects are normally created using the :py:func:`coldp.COLDP.start_name_bundle` method. 

At minimum a NameBundle is initialised with a dictionary holding a set of COLDP name record values. The scientific name represented by this dictionary should be the accepted name for a species or other taxon. Synonyms may then be added to the NameBundle. 

Once all names are included, the NameBundle can be added to the COLDP object via the :py:func:`coldp.COLDP.add_names` method. This adds name, taxon and synonym records for the set of names supplied. COLDP options may expand the set of added synonyms to include variant formats or may trigger the addition of one or taxa that are implicit in the accepted name. 

The COLDP object will automatically manage record identifiers and the basionymID for any name records that are combinations of another name in the set.

Installation
------------
Install the latest version from PyPI.

.. code-block:: console
    
    pip install py-coldp

Classes
-------
.. toctree::
   :maxdepth: 2

   coldp
   name-bundle
   identifier-policy

Usage
-----
.. toctree::
   :maxdepth: 2

   usage

Index and search 
~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`search`

coldp.COLDP
===========

.. currentmodule:: coldp

Class
~~~~~

.. autoclass:: COLDP
    :no-index:

Methods
~~~~~~~

Control behaviour
^^^^^^^^^^^^^^^^^
.. automethod:: coldp.COLDP.set_options
.. automethod:: coldp.COLDP.set_default_taxon_record
.. automethod:: coldp.COLDP.set_context

Add or modify records
^^^^^^^^^^^^^^^^^^^^^
.. automethod:: coldp.COLDP.add_names
.. automethod:: coldp.COLDP.modify_taxon
.. automethod:: coldp.COLDP.add_references
.. automethod:: coldp.COLDP.add_type_material
.. automethod:: coldp.COLDP.add_distribution
.. automethod:: coldp.COLDP.add_species_interaction
.. automethod:: coldp.COLDP.add_name_relation

Find or get records
^^^^^^^^^^^^^^^^^^^
.. automethod:: coldp.COLDP.get_reference

Tidy package
^^^^^^^^^^^^
.. automethod:: coldp.COLDP.fix_basionyms
.. automethod:: coldp.COLDP.fix_classification
.. automethod:: coldp.COLDP.sort_taxa
.. automethod:: coldp.COLDP.sort_names
.. automethod:: coldp.COLDP.reset_ids

Utilities
^^^^^^^^^
.. automethod:: coldp.COLDP.get_text_tree

Access to DataFrames
^^^^^^^^^^^^^^^^^^^^
.. automethod:: coldp.COLDP.table_by_name

Constants
~~~~~~~~~

.. autodata:: csv_extensions
    :no-value:
.. autodata:: id_mappings
    :no-value:
.. autodata:: name_from_nameusage
    :no-value:
.. autodata:: taxon_from_nameusage
    :no-value:
.. autodata:: synonym_from_nameusage
    :no-value:

Internal methods
~~~~~~~~~~~~~~~~

.. automethod:: coldp.COLDP.initialise_dataframe
.. automethod:: coldp.COLDP.extract_table

.. automethod:: coldp.COLDP.find_reference

.. automethod:: coldp.COLDP.fix_classification_recursive
.. automethod:: coldp.COLDP.sort_taxa_recursive

Index and search 
~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`search`
   
   
   
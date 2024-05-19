coldp.COLDP
===========

.. currentmodule:: coldp

Class
~~~~~

.. automethod:: coldp.COLDP.__init__

Methods
~~~~~~~

Control behaviour
^^^^^^^^^^^^^^^^^
.. automethod:: coldp.COLDP.set_options
.. automethod:: coldp.COLDP.set_default_taxon_record
.. automethod:: coldp.COLDP.set_context

Add or modify records
^^^^^^^^^^^^^^^^^^^^^
.. automethod:: coldp.COLDP.start_name_bundle
.. automethod:: coldp.COLDP.add_names
.. automethod:: coldp.COLDP.add_references
.. automethod:: coldp.COLDP.add_type_material
.. automethod:: coldp.COLDP.add_distribution
.. automethod:: coldp.COLDP.add_species_interaction
.. automethod:: coldp.COLDP.add_name_relation
.. automethod:: coldp.COLDP.add_synonym

.. automethod:: coldp.COLDP.modify_name
.. automethod:: coldp.COLDP.modify_taxon

Save
^^^^
.. automethod:: coldp.COLDP.save

Find or get records
^^^^^^^^^^^^^^^^^^^
.. automethod:: coldp.COLDP.find_taxon
.. automethod:: coldp.COLDP.find_name_record
.. automethod:: coldp.COLDP.find_names
.. automethod:: coldp.COLDP.find_name
.. automethod:: coldp.COLDP.find_reference
.. automethod:: coldp.COLDP.find_distribution
.. automethod:: coldp.COLDP.find_species_interaction
.. automethod:: coldp.COLDP.find_type_material
.. automethod:: coldp.COLDP.get_name
.. automethod:: coldp.COLDP.get_reference
.. automethod:: coldp.COLDP.get_taxon
.. automethod:: coldp.COLDP.get_synonyms
.. automethod:: coldp.COLDP.get_synonymy
.. automethod:: coldp.COLDP.get_children

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
.. automethod:: coldp.COLDP.get_available_column_headings
.. automethod:: coldp.COLDP.get_identifier_policy

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

.. automethod:: coldp.COLDP.insert_taxon
.. automethod:: coldp.COLDP.insert_synonym

.. automethod:: coldp.COLDP.fix_classification_recursive
.. automethod:: coldp.COLDP.sort_taxa_recursive

.. automethod:: coldp.COLDP.prepare_bundle
.. automethod:: coldp.COLDP.validate_record
.. automethod:: coldp.COLDP.identify_name
.. automethod:: coldp.COLDP.same_basionym
.. automethod:: coldp.COLDP.remove_gender
.. automethod:: coldp.COLDP.get_original_authorship
.. automethod:: coldp.COLDP.epithet_and_authorship_match
.. automethod:: coldp.COLDP.set_basionymid
.. automethod:: coldp.COLDP.fix_basionymid
.. automethod:: coldp.COLDP.construct_species_rank_name
.. automethod:: coldp.COLDP.construct_authorship
.. automethod:: coldp.COLDP.is_species_group
.. automethod:: coldp.COLDP.is_infrasubspecific
.. automethod:: coldp.COLDP.issue


Index and search 
~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`search`
   
   
   
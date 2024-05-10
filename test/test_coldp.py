#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Created By  : Donald Hobern, dhobern@gmail.com
# ---------------------------------------------------------------------------
""" 
Test suite for coldp.py classes
"""

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from coldp import COLDP, NameBundle
import pytest


def test_init():
    c = COLDP("Tonzidae", "input")
    print(c.names)
    t = c.find_taxon("Tonzidae", None, "family")
    print(t)
    print(c.get_text_tree(t["ID"]))
    assert len(c.names) == 5

import sys
import os
import pandas as pd
from coldp import COLDP

if __name__ == "__main__":
    if len(sys.argv) < 2:
        f"Usage: python {sys.argv[0]} coldp_folderpath"
        exit(0)

    coldp_folderpath = sys.argv[1]

    if not (os.path.exists(coldp_folderpath) and os.path.isdir(coldp_folderpath)):
        f"Error: {coldp_folderpath} is not a valid folder"
        exit(1)

    location, folder = os.path.split(coldp_folderpath)
    coldp = COLDP(location, folder, issues_to_stdout=True)

    if coldp is None:
        f"Error: Could not open {coldp_folderpath} as COLDP"
        exit(1)

    for n in coldp.names:
        match n:
            case { "genus": "Barea" }:
                f"{n}"

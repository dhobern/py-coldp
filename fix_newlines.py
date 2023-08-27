import csv
import os
import sys

if len(sys.argv) > 1:
    folder = sys.argv[1]
    for f in os.listdir(folder):
        filename = os.path.join(folder, f)
        if filename.endswith(".csv"):
            new_file = filename.replace(".csv", ".tsv")
            print(f"Processing {filename}")
            with open(filename, encoding="utf8", newline="") as infile, open(
                new_file, "w", encoding="utf8", newline=""
            ) as outfile:
                header = infile.readline()
                outfile.write(header)
                tabs = header.count("\t")
                print(f"{tabs + 1} columns")
                fixes = 0
                current = None
                for row in infile:
                    if current is not None:
                        row = current + " " + row
                        current = None
                    if row.count("\t") < tabs:
                        fixes += 1
                        current = row.replace("\n", "")
                    else:
                        outfile.write(row)
                print(f"{fixes} fixed")

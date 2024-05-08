import os
import sys
import re
from datetime import datetime

default_version = "0.1.000"
version_file = "VERSION"
version_pattern = re.compile(r"^.*\.[0-9]+$")
year_month_pattern = re.compile(r"^20[2-9][0-9]\.1?[0-9]$")


def get_version():
    if os.path.exists(version_file):
        with open(version_file) as f:
            version = f.readline().strip()
            if not version_pattern.match(version):
                version = default_version
    return version


def new_version():
    version = get_version()
    new_version = None

    parts = version.split(".")

    base = ".".join(parts[0 : len(parts) - 1])
    if year_month_pattern.match(base):
        new_base = datetime.today().strftime("%Y.%m").replace(".0", ".")
        if new_base != base:
            new_version = f"{new_base}.1"

    if new_version is None:
        number = parts[-1]
        length = len(number)
        new_version = f"{base}.{int(number)+1:0>{length}}"

    with open(version_file, "w") as f:
        f.write(new_version)

    return new_version


if __name__ == "__main__":
    if len(sys.argv) > 1:
        match sys.argv[1].lower():
            case "increment":
                new_version()
                exit(0)
            case "update":
                if len(sys.argv) == 5:
                    file = sys.argv[2]
                    if os.path.exists(file) and os.path.isfile(file):
                        version = get_version()
                        pattern = re.compile(sys.argv[3])
                        replacement = sys.argv[4].format(version=version)

                    with open(file) as f:
                        content = f.readlines()

                    with open(file, "w") as f:
                        for row in content:
                            f.write(pattern.sub(replacement, row))

                    exit(0)

print(f"python {sys.argv[0]} generate | update <file_name> <pattern> < replacement>\n")
print("  increment : Increment current version number in VERSION file")
print(
    "             If current version has form YYYY.MM.nnn, and the current date is not in YYYY.MM, the new version uses the current year and month and 001 as the final component."
)
print(
    "             Otherwise, the final component is incremented with a minimum of three digits and leading zeroes."
)
print(
    "  update : Replace supplied regex <pattern> in <filename> with the f-string <replacement>"
)
print(
    "             <replacement> should contain '{version}' to indicate where version number should appear"
)

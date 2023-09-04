# -*- coding: utf-8 -*-
from pathlib import Path
from zipfile import ZipFile
from majordome import download_file
import shutil

CASADI_VERSION = "3.6.1"
IPOPT_VERSION = "3.14.12"

HERE = Path(__file__).resolve().parent


def download_dependencies(deps):
    """ Manage download of all required files """
    # Construct path to releases page in GitHub.
    base_path_casadi = "https://github.com/casadi/casadi/releases/download"
    version_path_casadi = f"{base_path_casadi}/{CASADI_VERSION}"

    # List of required files to download.
    requires = {
        version_path_casadi : [
            f"casadi-cheatsheet_python-v{CASADI_VERSION}.pdf",
            f"casadi-users_guide-v{CASADI_VERSION}.pdf",
            f"casadi-example_pack-v{CASADI_VERSION}.zip",
            # f"casadi-source-v{CASADI_VERSION}.zip",
        ]
    }

    for pkg, files in requires.items():
        for fname in files:
            url = f"{pkg}/{fname}"
            saveas = deps / fname

            if saveas.exists():
                print(f"Skipping {saveas}, file exists")
                continue

            print(f"Downloading {fname}")
            download_file(url, saveas)


def extract_packages(deps):
    """ Extract zip containers inplace. """
    for zipname in deps.glob('*.zip'):
        with ZipFile(zipname) as zf:
            dest = deps / zipname.stem
            if dest.exists():
                print(f"Skipping {dest}, file exists")
                continue

            print(f"Extracting {zipname.stem}")
            zf.extractall(dest)


if __name__ == "__main__":
    # Ensure dependencies directory exists.
    deps = HERE / "deps"
    deps.mkdir(exist_ok=True)

    download_dependencies(deps)
    extract_packages(deps)

    dest = deps / f"casadi-source-v{CASADI_VERSION}/build.bat"
    shutil.copy("conf/build.bat", dest)

    print("Before running build manually flatten Ipopt directory")

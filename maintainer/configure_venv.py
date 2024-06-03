#pyMBE : The python-based Molecule Builder for ESPResSo.
#Copyright (C) 2023-2024 pyMBE-dev team
#
#This file is part of pyMBE.
#
#pyMBE is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#pyMBE is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import argparse
import sysconfig
try:
    import espressomd  # pylint: disable=unused-import
    espressomd_found = True
except ModuleNotFoundError:
    espressomd_found = False


def make_pth(name, path):
    if not os.path.isdir(path):
        raise ValueError(f"Folder '{path}' doesn't exist")
    site_packages = sysconfig.get_path("platlib")
    with open(os.path.join(site_packages, f"{name}.pth"), "w") as f:
        f.write(os.path.realpath(path))


parser = argparse.ArgumentParser(description="Configure pyBME and ESPResSo module paths")
parser.add_argument("--espresso_path", type=str, required=not espressomd_found,
                    help="Path to the ESPResSo build folder")
args = parser.parse_args()

if not os.environ.get("VIRTUAL_ENV"):
    raise RuntimeError("This script should be run in a virtual environment")

if not espressomd_found:
    make_pth("espresso", os.path.join(args.espresso_path, "src", "python"))
make_pth("pymbe", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

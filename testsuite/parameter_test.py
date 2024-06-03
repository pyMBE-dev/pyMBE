#pyMBE : The python-based Molecule Builder for ESPResSo.
#Copyright (C) 2024 pyMBE-dev team
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

import pathlib
import pyMBE

pmb = pyMBE.pymbe_library(SEED=42)

print("*** Check that the different pKa sets are correctly formatted ***")

pymbe_root = pathlib.Path(pyMBE.__file__).parent
pka_root = pymbe_root / "parameters" / "pka_sets"

for path in pka_root.glob("*.json"):
    print(f"Checking {path.stem}")
    path_to_pka = pmb.get_resource(path.relative_to(pymbe_root).as_posix())
    assert pathlib.Path(path_to_pka) == path
    pmb.load_pka_set(path_to_pka,
                     verbose=False)

print("*** Test passed ***")

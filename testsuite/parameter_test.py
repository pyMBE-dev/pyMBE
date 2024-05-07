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

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

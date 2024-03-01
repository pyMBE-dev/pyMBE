import os
import argparse
import sysconfig


def make_pth(name, path):
    if not os.path.isdir(path):
        raise ValueError(f"Folder '{path}' doesn't exist")
    site_packages = sysconfig.get_path("platlib")
    with open(os.path.join(site_packages, f"{name}.pth"), "w") as f:
        f.write(os.path.realpath(path))


parser = argparse.ArgumentParser(description="Configure pyBME and ESPResSo module paths")
parser.add_argument("espresso_build_path", type=str, help="Path to the ESPResSo build folder")
args = parser.parse_args()

if not os.environ.get("VIRTUAL_ENV"):
    raise RuntimeError("This script should be run in a virtual environment")

make_pth("pymbe", os.path.dirname(os.path.dirname(__file__)))
make_pth("espresso", os.path.join(args.espresso_build_path, "src", "python"))

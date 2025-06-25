# 
# Copyright (C) 2023-2024 pyMBE-dev team
#
# This file is part of pyMBE.
#
# pyMBE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyMBE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import stat
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


parser = argparse.ArgumentParser(description="Configure pyMBE and ESPResSo module paths")
parser.add_argument("--espresso_path", type=str, required=not espressomd_found,
                    help="Path to the ESPResSo build folder")
args = parser.parse_args()

if not os.environ.get("VIRTUAL_ENV"):
    raise RuntimeError("This script should be run in a virtual environment")

if not espressomd_found:
    make_pth("espresso", os.path.join(args.espresso_path, "src", "python"))
make_pth("pymbe", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def make_venv_set_omp(name, value):
    return f"""
_OLD_VIRTUAL_{name}="${name}"
if test -z "${name}"; then
  export {name}={value}
fi"""


def make_venv_unset_env(name):
    return f"""    if [ -n "${{_OLD_VIRTUAL_{name}:-}}" ] ; then
        {name}="${{_OLD_VIRTUAL_{name}:-}}"
        export {name}
        unset _OLD_VIRTUAL_{name}
    fi"""


if os.environ.get("VIRTUAL_ENV"):
    venv_activate = os.path.join(os.environ.get("VIRTUAL_ENV"), "bin", "activate")
    if os.path.isfile(venv_activate):
        original_access_mode = os.stat(venv_activate).st_mode
        os.chmod(venv_activate, original_access_mode | stat.S_IWUSR)
        try:
            with open(venv_activate, "r+") as f:
                content = f.read()
                token = "\nexport PATH"
                assert token in content
                content = content.replace(token, f"{token}\n{make_venv_set_omp('OMP_PROC_BIND', 'false')}\n{make_venv_set_omp('OMP_NUM_THREADS', '1')}")
                token = "unset _OLD_VIRTUAL_PATH\n    fi"
                assert token in content
                content = content.replace(token, f"{token}\n{make_venv_unset_env('OMP_PROC_BIND')}\n{make_venv_unset_env('OMP_NUM_THREADS')}""")
                f.seek(0)
                f.truncate()
                f.write(content)
        except:
            raise
        finally:
            os.chmod(venv_activate, original_access_mode)

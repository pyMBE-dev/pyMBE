#
# Copyright (C) 2024 pyMBE-dev team
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

# Create an instance of pyMBE library
import pyMBE
pmb = pyMBE.pymbe_library(seed=42)

print ('*** Unit tests: check that parse_sequence_from_file works correctly ***')

test_sequence = " , 'AAAAA','B ',' C', 'D''D', 'EEEEEE'"
reference_parsed_sequence = ["AAAAA", "B", "C", "DD"]
assert pmb.parse_sequence_from_file(test_sequence) == reference_parsed_sequence

print("*** Unit test passed***")

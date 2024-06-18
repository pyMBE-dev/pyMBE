# Please contribute to pyMBE!

## Bug reporting
If you discover any strange feature or unexpected behaviour of the library, please report it by opening an issue in this GitHub repository or contacting Dr. Pablo M. Blanco (pablb@ntnu.no) or any other member of the pyMBE development team.
Once a ticket is open in GitHub, we will work to fix the issue as soon as possible.

## New contributors
New developers are welcome to contribute to extend the functionalities of the pyMBE library.
To contribute to the pyMBE library, first one needs to be added as a member of this GitHub repository.
If you want to contribute to the development of the pyMBE library, please contact Dr. Pablo M. Blanco (pablb@ntnu.no)

## Rules to contribute to pyMBE
Create a fork of the repository and submit your contribution in the form of a pull request.
Give your branch a short and meaningful name, and start your work from the most recent commit on the main branch.
Any new version of the code must reproduce all the data stored in `testsuite/data`.
Before pushing your code, run `make tests` to execute the full test suite,
`make pylint` to check for code issues, and `make docs` to confirm sphinx can build the user guide.
When rapidly prototyping code, run `make unit_tests` periodically to check for regressions.
All new code will be reviewed by at least one member of the pyMBE development team before being merged into the main branch to ensure that a functional version of the code is always available.
Class methods are sorted in alphabetical order.
We follow semantic versioning and keep a changelog.

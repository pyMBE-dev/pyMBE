# Please contribute to pyMBE!

## Bug reporting
If you discover any strange feature or unexpected behaviour of the module, please report it by opening an issue in this GitHub repository or contacting Dr. Pablo M. Blanco (pablb@ntnu.no) or any other member of the pyMBE development team.
Once a ticket is open in GitHub, we will work to fix the issue as soon as possible.

## New contributors
New developers are welcomed to contribute to extend the functionalities of the pyMBE module. 
To contribute to the pyMBE module, first one needs to be added as a member of this GitHub repository.
If you want to contribute to the development of the pyMBE module, please contact Dr. Pablo M. Blanco (pablb@ntnu.no).

## Authorship system
We acknowledge as valid contributions to the module:
- New methods and modules that extend the functionalies of the module.
- Bugfixes and code enhancements that improve the quality of the module.
- Development of tutorials and other documentation of the module.
- Scripts for samples on how to use pyMBE, testing of the library or as handy tools to extend the features fo the library.
- Ideas on the software design and carpentry. This includes ideas on how to improve the code quality, ideas on how to improve the existing features of the module, ideas for adding new features to the module and other conceptual contributions that have an impact on the positive development of pyMBE.

In the pyMBE community, we use a tier system to select how to reward our contributors: 

### Tier-0 
- Contributors of big features of pyMBE, including stand-alone modules or a set of several methods.
- Admins of pyMBE.
- Contributors whose ideas had a substancial impact on the carpentry of the library.

### Tier-1
- Occasional developers of methods for pyMBE.
- Contributors of bugfixes and small code enhancements.
- Contributors of single-purpose scripts and tutorials.

All tier contributors will be acknowledged as authors of the module in `AUTHORS.md`. 
Additionally, Tier-0 contributors will be acknowledged in the `CITATION.cff`.
New Tier-0 contributors will be included in the `CITATION.cff` file after the corresponding release of a new version of the software.

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

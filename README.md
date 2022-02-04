# CombiFF

![combiff](https://user-images.githubusercontent.com/13115540/139452181-34123bf6-68ad-45d0-a0f7-662f15ad29fb.png)

## Introduction

## Usage

### Installation

For installation instructions please refer to [INSTALL.md](https://github.com/csms-ethz/CombiFF/blob/main/INSTALL.md)

### Getting Started

#### Scripts

There are several scripts that can help you with the use of CombiFF. The default settings are defined in `use/global_settings.py`. This file should *not* be modified to avoid merge conflicts. Instead, you can define your own local settings in `use/user_settings.py` which is ignored by git. In this file, you can re-define any of the variables defined in the global settings file to override the default choices.

 * `scr/clean.py`: clean up of temporary files
 * `scr/compile.py`: compile `enu`, `tbl`, and `cnv`.
 * `scr/compile_documentation.py`: generate PDFs of the documentation from the tex files in doc
 * `run_enu.py`: run `enu` with the specifications of `use/global_settings.py` and/or `use/user_settings.py`
 * `run_tbl.py`: run `tbl` with the specifications of `use/global_settings.py` and/or `use/user_settings.py`
 * `draw_isomers.py`: draw enumerated isomers with RDKit

Known variables for the scripts are:

 * `families_to_enumerate`: list of strings, containing the code of the families that should be enumerated by `enu`. If it contains the element `all`, all families will be enumerated. You can use regex such as `'01..'` or `*1*`.
 * `families_to_build_topologies`: list of strings, containing the code of the families for which topologies should be generated by `tbl`. If it contains the element `all`, topologies will be generated for all families. You can use regex such as `'01..'` or `*1*`.
 * `family_exceptions`: list of strings, containing the code of the families that should be ignored by the scripts, even if `all` is specified for `families_to_enumerate` or `families_to_build_topologies`. You can use regex such as `'01..'` or `*1*`.
 * `force_update`: boolean, specifying whether output files should be regenerated even if input files did not change since the last execution of `enu` or `tbl`
 * `united_atoms`: boolean, specifying whether all atom or united atom convention should be used for aliphatic CHn groups
 * `third_neighbor_exclusions`: boolean, specifying whether 3rd neighbors should also be excluded
 * `unique_torsionals`: boolean, specifying whether only one torsional per central bond should be reported in a topology
 * `stereo`: boolean, specifying whether stereoisomers should also be enumerated by `enu`
 * `fragments_to_use`: list of strings, containing the fragments that should be used by `tbl`, in order of decreasing priority. If it is set to `None`, all fragments will be used, and the priority is determined according to priority rules (see documentation of `tbl`).

#### First Steps

When you familiarize with the setup of CombiFF, a good starting place is the global settings file `use/global_settings.py`. Here you will see the definition of the default input and output directories, as well as the default settings for `enu` and `tbl`. By default, the family with the code `intro` will be used for both `enu` and `tbl`. It is defined in `use/input_files/family_definitions/family_definitions_intro.xml`. You can create a list of all conformational isomers of the `intro` family by executing the script `scr/run_enu.py`. This will create the xml file `use/output_files/family_isomer_enumeration_intro.xml`. Afterwards, you can generate a GROMOS mtb file by executing the script `scr/run_tbl.py`. This will create intermediate files in `use/output_files/molecule_decompositions` (specifying which fragments are linked together to create the molecules) and in `use/output_files/molecules_with_macros` (topology in xml format), and it will create an mtb file in `use/output_files/mtb` (note: currently still with macro parameters).

Next, you could have a look at the different families defined in `use/input_files/family_definitions` and play around with `enu` and `tbl` by adapting your settings in `use/user_settings.py` to e.g. choose different families, or enumerate stereo isomers. You can have a look at the `README` files within the different directories of `use/input_files` if you want to define your own families, fragments, etc. Ideally you shouldn't modify the existing input files but rather create your own. Just follow the format outlined in the corresponding `README` and `.dtd` files.

You can of course also work with the differen programs without the python helper scripts. Simply start by executing them with the `--help` option to see a list of options, or have a look at the documentation.

### Enu

For documentation of `enu`, please refer to `doc/enu.pdf` (coming soon). Please make sure to compile the documentation first by executing the script `scr/compile_documentation.py`

### Tbl

For documentation of `tbl`, please refer to `doc/tbl.pdf` (coming soon). Please make sure to compile the documentation first by executing the script `scr/compile_documentation.py`

### Cnv

For documentation of `cnv`, please refer to `doc/cnv.pdf` (coming soon). Please make sure to compile the documentation first by executing the script `scr/compile_documentation.py`

## References

Marina P. Oliveira, Maurice Andrey, Salomé R. Rieder, Leyla Kern, David F. Hahn, Sereina Riniker, Bruno A. C. Horta, and Philippe H. Hünenberger - Systematic Optimization of a Fragment-Based Force Field against Experimental Pure-Liquid Properties Considering Large Compound Families: Application to Saturated Haloalkanes. J. Chem. Theory Comput. 2020, 16, 12, 7525–7555. https://doi.org/10.1021/acs.jctc.0c00683

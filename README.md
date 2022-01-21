# CombiFF

![combiff](https://user-images.githubusercontent.com/13115540/139452181-34123bf6-68ad-45d0-a0f7-662f15ad29fb.png)

## Introduction

## Usage

For installation instructions please refer to INSTALL.md

### Scripts

Default settings are defined in use/global_settings.py. This file should *not* be modified to avoid merge conflicts. Instead, you can define your own local settings in use/user_settings.py which is ignored by git.

There are several scripts that can help you with the use of CombiFF:

 * scr/clean.py: clean up of temporary files
 * scr/compile.py: compile enu, tbl, and cnv
 * scr/compile_documentation.py: generate PDFs of the documentation from the tex files in doc
 * draw_isomers.py: draw enumerated isomers with RDKit
 * run_enu.py: run enu with the specifications of use/global_settings.py and/or use/user_settings.py
 * run_tbl.py: run tbl with the specifications of use/global_settings.py and/or use/user_settings.py

### Enu

### Tbl

### Cnv

## Examples

## References

Marina P. Oliveira, Maurice Andrey, Salomé R. Rieder, Leyla Kern, David F. Hahn, Sereina Riniker, Bruno A. C. Horta, and Philippe H. Hünenberger - Systematic Optimization of a Fragment-Based Force Field against Experimental Pure-Liquid Properties Considering Large Compound Families: Application to Saturated Haloalkanes. J. Chem. Theory Comput. 2020, 16, 12, 7525–7555. https://doi.org/10.1021/acs.jctc.0c00683

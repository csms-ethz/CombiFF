import os
import re

#
# Important:
# you should not change the default choices in this file to avoid conflicts with different versions of other users on git
# if you'd like to make a different choice, please create a file 'user_settings.py' and modify your parameters in there
#

root_dir = os.path.abspath(os.path.dirname(__file__) + "/..")

def build_path(directory, subdirectory):
  return directory + '/' + subdirectory

#
# directory definitions
#

doc_dir = build_path(root_dir, 'doc')
use_dir = build_path(root_dir, 'use')
input_dir = build_path(use_dir, 'input_files')
output_dir = build_path(use_dir, 'output_files')
build_dir = build_path(root_dir, 'build')
bin_dir = build_path(root_dir, 'bin')

enu_executable = build_path(bin_dir, 'enu')
tbl_executable = build_path(bin_dir, 'tbl')
cnv_executable = build_path(bin_dir, 'cnv')

aliases_dir = build_path(input_dir, 'aliases')
atomtypes_dir = build_path(input_dir, 'atomtypes')
family_definitions_dir = build_path(input_dir, 'family_definitions')
fragments_dir = build_path(input_dir, 'fragments')
pseudoatoms_dir = build_path(input_dir, 'pseudoatoms')
substructures_dir = build_path(input_dir, 'substructures')

family_isomer_enumerations_dir = build_path(output_dir, 'family_isomer_enumerations')
family_isomer_enumerations_stereo_dir = build_path(output_dir, 'family_isomer_enumerations_stereo')
molecule_decompositions_dir = build_path(output_dir, 'molecule_decompositions')
molecules_with_macros_dir = build_path(output_dir, 'molecules_with_macros')
mtb_dir = build_path(output_dir, 'mtb')
img_dir = build_path(output_dir, 'img')

#
# default choices
#

family_exceptions=['test*']
families_to_enumerate = ['0100']
families_to_build_topologies = ['0100']

force_update = True
united_atoms = True
third_neighbor_exclusions = True
unique_torsionals = True
stereo = False
fragments_to_use = None #leave empty to use all fragments

#
# overwrite default choices with user settings
#

if(os.path.exists(os.path.abspath(os.path.dirname(__file__)) + '/user_settings.py')):
  from user_settings import *

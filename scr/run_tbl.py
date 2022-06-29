import sys
import subprocess
import glob
from pathlib import Path
import filecmp
import xml.etree.ElementTree as ET
from helper_functions import generate_backup_path
from helper_functions import elements_equal
import shutil
import os

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(script_dir + "/../use")

from global_settings import *

os.makedirs(molecule_decompositions_dir, exist_ok=True)
os.makedirs(molecules_with_macros_dir, exist_ok=True) 
os.makedirs(mtb_dir, exist_ok=True) 

enumerations = glob.glob(family_isomer_enumerations_dir + "/*.xml")
enumerations.sort()
enumerations_to_update = []

for enumeration in enumerations:
  tree_enumeration = ET.parse(enumeration)
  root_enumeration = tree_enumeration.getroot()
  enumeration_definitions = root_enumeration.findall('isomer_lists')

  if not Path(generate_backup_path(enumeration)).exists() or force_update:
    enumerations_to_update.append((root_enumeration.attrib['family_code'], enumeration))

  else:
    tree_backup = ET.parse(generate_backup_path(enumeration))
    root_backup = tree_backup.getroot()

    backup_enumeration_definitions = root_backup.findall('isomer_lists')

    for enum, backup_enum in zip(enumeration_definitions, backup_enumeration_definitions):
      if not elements_equal(enum, backup_enum):
        enumerations_to_update.append((root_enumeration.attrib['family_code'], enumeration))
      elif not len(glob.glob(build_path(molecule_decompositions_dir, '*' + root_enumeration.attrib['family_code'] + '.xml'))):
        enumerations_to_update.append((root_enumeration.attrib['family_code'], enumeration))


enumerations_to_update = dict(enumerations_to_update)
print("\nenumerations to update: ", enumerations_to_update.keys())

family_exceptions_re=re.compile(r'\b(?:%s)\b' % '|'.join(family_exceptions))
enumerations_to_update = [(fam, fil) for fam, fil in enumerations_to_update.items() if not family_exceptions_re.match(fam)]
enumerations_to_update = dict(enumerations_to_update)

print("\nremoving exceptions: ", family_exceptions)
print("enumerations to update: ", enumerations_to_update.keys())

families_to_build_topologies_re=re.compile(r'\b(?:%s)\b' % '|'.join(families_to_build_topologies))
enumerations_to_update = [(fam, fil) for fam, fil in enumerations_to_update.items() if families_to_build_topologies_re.match(fam) or 'all' in families_to_build_topologies]
enumerations_to_update = dict(enumerations_to_update)

print("\nchecking for explicit family list: ", families_to_build_topologies)
print("enumerations to update: ", enumerations_to_update.keys())

fragment_files = glob.glob(fragments_dir + "/*.xml")
if (not united_atoms):
    fragment_files = [ff for ff in fragment_files if "united" not in ff]

replacement_files = glob.glob(replacement_dir + "/*.xml")

if(not os.path.exists(tbl_executable)):
  print("Warning: tbl executable not found -> compiling")
  exec(open(script_dir + '/compile.py').read())

if len(enumerations_to_update):
  for family, file in enumerations_to_update.items():
    print(family)
    
    command = [tbl_executable, "-fragments"]
    command.extend(fragment_files)
    command.append("-replacements")
    command.extend(replacement_files)
    command.extend(["-output_directory_decompositions", molecule_decompositions_dir])
    command.extend(["-output_directory_molecule_macros", molecules_with_macros_dir])
    command.extend(["-output_directory_mtb", mtb_dir])
    command.append("-families")
    command.append(file)
    command.append("-atomtypes")
    command.extend(glob.glob(atomtypes_dir + "/*.xml"))
    if(united_atoms):
      command.append("-united_atom")
      
    if(fragments_to_use is not None and len(fragments_to_use)):
      command.append("-fragments_to_use")
      command.extend(fragments_to_use)
      
    if(third_neighbor_exclusions):
      command.append("-third_neighbor_exclusions")
      
    if(unique_torsionals):
      command.append("-unique_torsionals")

    print(command)
    try:
      pr = subprocess.run(command, stderr=subprocess.PIPE)
      pr.check_returncode()
      
    except subprocess.CalledProcessError as e:
      print ( "Error:\nreturn code: ", e.returncode, "\nOutput: ", e.stderr.decode("utf-8") )
      raise
    
    shutil.copyfile(file, generate_backup_path(file))

else:
  print("everything up to date")

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


print(family_isomer_enumerations_dir)
os.makedirs(family_isomer_enumerations_dir, exist_ok=True) 
os.makedirs(family_isomer_enumerations_stereo_dir, exist_ok=True)

families = glob.glob(family_definitions_dir + "/*.xml")

families.sort()
families_to_update = []

for family in families:
  tree_family = ET.parse(family)
  root_family = tree_family.getroot()
  family_definitions = root_family.findall('family_definition')

  if not Path(generate_backup_path(family)).exists() or force_update:
    families_to_update.extend([(fam.attrib['code'], family) for fam in family_definitions])

  else:
    tree_backup = ET.parse(generate_backup_path(family))
    root_backup = tree_backup.getroot()

    backup_family_definitions = root_backup.findall('family_definition')

    for fam, fam_backup in zip(family_definitions, backup_family_definitions):
      if not elements_equal(fam, fam_backup):
        families_to_update.append((fam.attrib['code'], family))
      elif not len(glob.glob(build_path(family_isomer_enumerations_dir, '*' + fam.attrib['code'] + '.xml'))):
        families_to_update.append((fam.attrib['code'], family))

families_to_update = dict(families_to_update)
print("\nfamilies to update: ", families_to_update.keys())

family_exceptions_re=re.compile(r'\b(?:%s)\b' % '|'.join(family_exceptions))
families_to_update = [(fam, fil) for fam, fil in families_to_update.items() if not family_exceptions_re.match(fam)]
families_to_update = dict(families_to_update)

print("\nremoving exceptions: ", family_exceptions)
print("familes to update: ", families_to_update.keys())

families_to_enumerate_re=re.compile(r'\b(?:%s)\b' % '|'.join(families_to_enumerate))
families_to_update = [(fam, fil) for fam, fil in families_to_update.items() if families_to_enumerate_re.match(fam) or 'all' in families_to_enumerate]
families_to_update = dict(families_to_update)

print("\nchecking for explicit family list: ", families_to_enumerate)
print("families to update: ", families_to_update.keys())

if(not os.path.exists(enu_executable)):
  print("Warning: enu executable not found -> compiling")
  exec(open(script_dir + '/compile.py').read())

if len(sys.argv) > 1:
  command = sys.argv[1:len(sys.argv)]
  command.insert(0, enu_executable)
  subprocess.call(command)
  
elif len(families_to_update):
  for family, file in families_to_update.items():
    print(family)
    
    command = [enu_executable, "-family_files"]
    command.append(file)
    command.append("-families")
    command.append(family)
    if(not stereo):
        command.extend(["-output_directory", family_isomer_enumerations_dir])
    else:
    	command.extend(["-output_directory", family_isomer_enumerations_stereo_dir])
    	command.append("-stereo")
    command.append("-substructure_files")
    command.extend(glob.glob(substructures_dir + "/*.xml"))
    command.append("-element_alias_files")
    command.extend(glob.glob(aliases_dir + "/*.xml"))
    command.append("-pseudoatom_files")
    command.extend(glob.glob(pseudoatoms_dir + "/*.xml"))
    print(command)
    success = not subprocess.call(command)
    if(success):
      shutil.copyfile(file, generate_backup_path(file))
    else:
      print("ERROR")
      exit()

else:
  print("everything up to date")


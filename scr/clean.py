import os
import sys
import glob

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(script_dir + "/../use")

from global_settings import *

os.chdir(family_definitions_dir)
for fil in glob.glob(".*.xml*"):
  os.remove(fil)


os.chdir(family_isomer_enumerations_dir)
for fil in glob.glob(".*.xml*"):
  os.remove(fil)

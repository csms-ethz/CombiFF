import sys
import subprocess
import glob
import filecmp
import xml.etree.ElementTree as ET
from helper_functions import generate_backup_path
from helper_functions import elements_equal
import shutil
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('executable', default = None, nargs='?')

args = parser.parse_args()

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(script_dir + "/../use")

from global_settings import *

os.makedirs(build_dir, exist_ok=True)
os.chdir(build_dir)

subprocess.call(['cmake','..'])
  
if(args.executable is not None):
  subprocess.call(['make','-j4',args.executable])
else:
  subprocess.call(['make','-j4'])
  #copy script files
  for file in glob.glob(script_dir + "/*.py"):
    shutil.copy(file, bin_dir)
  

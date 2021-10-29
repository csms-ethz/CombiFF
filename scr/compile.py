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

os.makedirs(build_dir, exist_ok=True)
os.chdir(build_dir)

subprocess.call(['cmake','..'])
  
subprocess.call(['make','-j4'])

import os
import sys
import filecmp
import shutil
from pathlib import Path

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(script_dir + "/../use")

from global_settings import doc_dir

def recompile(cff_man, cff_man_bu):
  print("recompile")
  os.system('./jcompile')
  shutil.copyfile(cff_man, cff_man_bu)

os.chdir(doc_dir)

cff_man = Path('cff_man.tex')
cff_man_bu = Path('.cff_man.tex')

if not cff_man.exists():
  print("cff_man.tex not found")
  exit()

elif not cff_man_bu.exists() or not filecmp.cmp('cff_man.tex', '.cff_man.tex'):
  recompile(cff_man = cff_man, cff_man_bu = cff_man_bu)

else:
  print("cff_man is up to date")

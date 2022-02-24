import os
import sys
import filecmp
import shutil
import glob
import subprocess

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(script_dir + "/../use")

from global_settings import doc_dir

def recompile(cff_man, cff_man_bu):
  print("recompile")
  os.system('./jcompile')
  shutil.copyfile(cff_man, cff_man_bu)

os.chdir(doc_dir + "/tex")

for tex_file in glob.glob(doc_dir + "/tex/*.tex"):
  doc_file = os.path.basename(tex_file)
  doc_file_compiled = "." + doc_file
  
  pdf_file = glob.glob(doc_dir + "/" + "".join(doc_file.split('.')[:-1]) + ".pdf")
  print(pdf_file)
  print(doc_file)  

  if not len(pdf_file) or not os.path.isfile(doc_file_compiled) or not filecmp.cmp(doc_file, doc_file_compiled):
    if(len(pdf_file) == 1):
      os.remove(pdf_file[0])
    elif(len(pdf_file) > 1):
      print("ERROR: too many files found: ", pdf_file)
    subprocess.run(["pdflatex", doc_file])
    pdf_file = glob.glob("".join(doc_file.split('.')[:-1]) + ".pdf")[0]
    shutil.copyfile(doc_file, doc_file_compiled)
    shutil.move(pdf_file, doc_dir)

  else:
    print(doc_file + " is up to date")

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing
import xml.etree.ElementTree as ET
from pdfrw import PdfReader, PdfWriter
import os
import sys
import glob

script_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(script_dir + "/../use")

from global_settings import *

if(len(sys.argv) < 1):
  print("please provide input file")
  exit()

input_file = sys.argv[1]

input_file = open(input_file)

lines = input_file.readlines()
for index, line in enumerate(lines):
  smiles = [s for s in line.split(' ') if s != "" and s != "\n"]
  if(len(smiles) != 2):
    continue
  print(smiles)
  mol1 = Chem.MolFromSmiles(smiles[0])
  mol2 = Chem.MolFromSmiles(smiles[1])
  if(Chem.MolToSmiles(mol1) != Chem.MolToSmiles(mol2)):
    print(smiles, "not equal!!")

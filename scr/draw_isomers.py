from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
import xml.etree.ElementTree as ET
from pdfrw import PdfReader, PdfWriter
import os
import sys

DrawingOptions.atomLabelFontSize = 50
DrawingOptions.dotsPerAngstrom = 100
DrawingOptions.bondLineWidth = 2.0

if(len(sys.argv) < 1):
  print("please provide family isomer enumeration file")
  exit()

tree = ET.parse(sys.argv[1])
root = tree.getroot()

mol_list = []
id_list = []

num_iso = root.find('number_of_isomers').text
num_printed = 0

molsPerRow = 10
molsPerPage = int(molsPerRow*(molsPerRow)*2)
writer = PdfWriter("isomerlist.pdf")

for isomer_lists in root.findall('isomer_lists'):
  for isomer_list in isomer_lists.findall('isomer_list'):
    for isomer in isomer_list.findall('isomer'):
      if len(mol_list) < molsPerPage:
        if isomer.find('stereo_isomers') is None:
          #print(isomer.find('constitutional_SMILES').text)
          mol_list.append(Chem.MolFromSmiles(isomer.find('constitutional_SMILES').text))
          #id_list.append(isomer.get('isomer_id'))
          id_list.append(isomer.find('constitutional_SMILES').text)
        else:
          for stereo_isomers in isomer.findall('stereo_isomers'):
            for stereo_isomer in stereo_isomers.findall('stereo_isomer'):
              if len(mol_list) < molsPerPage:
                #print(stereo_isomer.find('stereo_SMILES').text)
                mol_list.append(Chem.MolFromSmiles(stereo_isomer.find('stereo_SMILES').text))
                #id_list.append(stereo_isomer.get('stereo_id'))
                id_list.append(stereo_isomer.find('stereo_SMILES').text)
              else:
                 img = Chem.Draw.MolsToGridImage(mol_list, molsPerRow=molsPerRow, subImgSize = (200,200), maxMols = molsPerPage, legends=id_list).convert('RGB')
                 img.save("tmp.pdf")
                 writer.addpages(PdfReader("tmp.pdf").pages)
                 mol_list.clear()
                 id_list.clear()
                 num_printed += molsPerPage
                 print(str(num_printed) + "/" + str(num_iso))
                
      else:
        img = Chem.Draw.MolsToGridImage(mol_list, molsPerRow=molsPerRow, subImgSize = (200,200), maxMols = molsPerPage, legends=id_list).convert('RGB')
        img.save("tmp.pdf")
        writer.addpages(PdfReader("tmp.pdf").pages)
        writer.write()
        mol_list.clear()
        id_list.clear()
        num_printed += molsPerPage
        print(str(num_printed) + "/" + str(num_iso))


for i in range(len(mol_list)):
  for j in range(i+1, len(mol_list)):
    m1 = mol_list[i]
    m2 = mol_list[j]
    s1 = Chem.CanonSmiles(Chem.rdmolfiles.MolToSmiles(m1))
    s2 = Chem.CanonSmiles(Chem.rdmolfiles.MolToSmiles(m2))
    print(str(i) + " " + str(j) + " " + str(s1 == s2))

img = Chem.Draw.MolsToGridImage(mol_list, molsPerRow=molsPerRow, subImgSize = (200,200), maxMols = molsPerPage, legends=id_list).convert('RGB')
img.save("tmp.pdf")
writer.addpages(PdfReader("tmp.pdf").pages)

writer.write()
os.remove("tmp.pdf")

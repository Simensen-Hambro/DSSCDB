#!/usr/bin/python

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


if len(sys.argv) != 4:
    sys.exit('Invalid format, format needs to be: SMILES output_file file_format')



smiles = sys.argv[1]
outfile = sys.argv[2]
format = sys.argv[3]
try:
    wrds = smiles.split(" ")
    if len(wrds) == 1:
        mol = Chem.MolFromSmiles(smiles)
    else:
        mol = Chem.MolFromSmiles(wrds[0])
    AllChem.Compute2DCoords(mol)
    image_url = "molecules/" + outfile + "." + format
    image_store_url = "media/" + image_url
    # image = Draw.MolsToImage([mol],subImgSize=(500,500))
    Draw.MolToFile(mol, image_store_url, size=(400, 400), type=format)
    sys.stdout.write('Image {} generated for SMILES: {}'.format(image_url, smiles))
except Exception as e:
    sys.stdout.write('Error message: {}, for SMILES: {}'.format(e, smiles))

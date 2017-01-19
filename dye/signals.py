from dye.models import Molecule
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.core.management import call_command
from uuid import uuid4
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

@receiver(post_save, sender=Molecule)
def generate_image(sender, instance, signal, created, **kwargs):
    if created:
        image_name = uuid4().hex[:10]
        smiles = instance.smiles
        format = 'png'
        try:
            wrds = smiles.split(" ")
            if len(wrds) == 1:
                mol = Chem.MolFromSmiles(smiles)
            else:
                mol = Chem.MolFromSmiles(wrds[0])
            AllChem.Compute2DCoords(mol)
            image_url = "molecules/" + image_name + "." + format
            image_store_url = "static/" + image_url
            Draw.MolToFile(mol, image_store_url, size=(500, 500), type=format)
            instance.image = image_url
            instance.save()
        except Exception as e:
            print('Failed to generate image, {}'.format(e))




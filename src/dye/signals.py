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
        call_command('generateimage', instance.pk, smiles, 'svg')




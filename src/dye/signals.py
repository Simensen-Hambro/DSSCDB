import logging
import os
from uuid import uuid4

from django.db.models.signals import post_save
from django.dispatch import receiver
from dye.models import Molecule
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from django.conf import settings
from .helpers import generate_coordinates_babel
import pybel


@receiver(post_save, sender=Molecule)
def generate_molecule_data(sender, instance, signal, created, **kwargs):
    """
    A signal to create extra data for molecules.
    2D image representation, 3D-coordinates (SDF) and fingerprint string for
    search.
    """
    logger = logging.getLogger('image')
    if created:

        # SVG file generation
        format = kwargs.get('format') or 'svg'

        image_name = uuid4().hex[:10]
        smiles = instance.smiles
        molecule = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(molecule)
        file_name = 'molecules/' + '{}.{}'.format(image_name, format)
        image_url = settings.MEDIA_ROOT + '/' +  file_name
        # 2D representation by RDKIT
        try:
            os.makedirs(os.path.dirname(image_url), exist_ok=True)
            Draw.MolToFile(molecule, image_url, size=(400, 400), type=format)
            instance.image = file_name
            logger.info('INFO: Image {} generated for SMILES: {}'.format(image_url, smiles))
        except Exception as e:
            logger.error('ERROR: Error message: {}, for SMILES: {}'.format(e, smiles))
        # SDF (3D-coordinates) generations
        sdf_string = generate_coordinates_babel(instance.smiles)
        if sdf_string == '':
            logger.error('ERROR: Failed to generate SDF for {}'.format(instance.smiles))

        instance.representation_3d = sdf_string

        pybel_object = pybel.readstring('smiles', instance.smiles)
        instance.fingerprint = str(pybel_object.calcfp(fptype='fp2'))

        instance.save()

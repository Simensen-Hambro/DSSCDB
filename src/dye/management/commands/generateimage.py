import sys
import subprocess
from django.conf import settings
from django.core.management.base import BaseCommand
from uuid import uuid4
from datetime import datetime
from dye.models import Molecule
# import the logging library
import logging


class Command(BaseCommand):
    args = 'smile outfile format'
    help = 'Create images of molecule from SMILE string'

    def add_arguments(self, parser):
        parser.add_argument('pk', type=int)
        parser.add_argument('smiles')
        parser.add_argument('format')


    def handle(self, *args, **options):
        smiles = options.get('smiles')
        format = options.get('format')
        pk = options.get('pk')
        generate_image(pk, smiles, format)


def generate_image(pk, smiles, format):
    logger = logging.getLogger('image')
    python_version = settings.PYTHON_2_ENV
    script_dir = settings.GENERATE_IMAGE_SCRIPT
    smiles = smiles
    out_file = uuid4().hex[:10]
    format = format
    process = subprocess.Popen([python_version, script_dir, smiles, out_file, format], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, )
    output, error = process.communicate()

    if output:
        output = output.decode('utf-8')

    if error:
        error = error.decode('utf-8')
    now = datetime.now()
    if error:
        logger.error('(ERROR {}): {}'.format(now, error))
    if output:
        logger.info('(INFO {}): {}'.format(now, output))
        molecule = Molecule.objects.get(pk=pk)
        molecule.image = 'molecules/' + out_file + '.' + format
        molecule.save()


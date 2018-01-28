from __future__ import unicode_literals

from concurrent import futures
from django.db import migrations
from django.db.models.signals import post_save
from ..helpers import generate_coordinates_babel



def reformat(app, schema_editor):
    try:
        Molecule = app.get_model('dye', 'Molecule')
    except LookupError:
        print("Old model is no longer installed")
        return
    post_save.disconnect(Molecule, sender=Molecule)

    molecules = Molecule.objects.all()
    with futures.ProcessPoolExecutor(max_workers=8) as executor:
        future_generate_coordinates = {executor.submit(generate_coordinates_babel, molecule.smiles): molecule for molecule
                                       in
                                       molecules}
        for future in futures.as_completed(future_generate_coordinates):
            molecule = future_generate_coordinates[future]
            data = future.result()
            molecule.representation_3d = data
            molecule.save()


class Migration(migrations.Migration):
    operations = [
        migrations.RunPython(reformat, migrations.RunPython.noop),
    ]
    dependencies = [
        ('dye', '0003_auto_20180128_1143'),
    ]

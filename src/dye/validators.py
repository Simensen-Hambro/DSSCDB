from django.core.exceptions import ValidationError
from pybel import readstring
from rdkit.Chem import MolFromInchi


def validate_smiles(data):
    molecule = None
    try:
        molecule = readstring('smiles', data)

    # Catch all, as Pybel doesn't have useful Python exceptions
    except:
        pass

    if not molecule:
        raise ValidationError('"{}" is not a valid SMILES'.format(data))


def validate_inchi(data):
    molecule = MolFromInchi(data)
    if not molecule:
        raise ValidationError('"{}" is not a valid INCHI'.format(data))


def not_negative(data):
    return data > 0


def between_0_and_1(data):
    return 0 < data < 1

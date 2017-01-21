from rdkit.Chem import MolFromInchi, MolFromSmiles


def validate_smiles(data):
    molecule = MolFromSmiles(data)
    if molecule:
        return True
    else:
        return False

def validate_inchi(data):
    molecule = MolFromInchi(data)
    if molecule:
        return True
    else:
        return False

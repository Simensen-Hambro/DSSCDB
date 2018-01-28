import re

import bibtexparser
import requests
from django.utils.timezone import datetime
import pybel


def get_DOI_metadata(doi):
    url = 'http://dx.doi.org/' + doi
    headers = {'accept': 'application/x-bibtex', 'style': 'bibtex'}
    response_bytes = requests.get(url, headers=headers).content
    response_string = response_bytes.decode("utf-8")
    bibtex_object = bibtexparser.loads(response_string)

    try:
        article = bibtex_object.entries[0]
    except IndexError:
        return None

    new_article_data = {
        'author': article.get('author'),
        'title': article.get('title'),
        'journal': article.get('journal'),
        'volume': article.get('volume'),
        'doi': article.get('doi'),
        'pages': article.get('pages'),
        'electronic_id': article.get('ID'),
        'issue_nr': article.get('number'),
        'keywords': article.get('keywords'),
        'year': datetime(year=int(article.get('year')), month=1, day=1),
    }

    return new_article_data


def to_decimal(string):
    # Shapes "string" into a number with a best-effort attempt
    if not isinstance(string, float):
        illegal_characters = re.search('([^0-9^.^,^-])', string)
        if illegal_characters:
            return string
        else:
            rep = re.compile('(\-?\d*\.?\d+)')

            result = rep.search(string)
            if result:
                return result
            else:
                return None
    else:
        if len(str(string)) >= 7:
            string = round(string, 6)
        return string


def locate_start_data(sheet):
    """
    Locate control tag start row
    """
    start_data = None
    for row_index in range(sheet.nrows):
        if "**%BEGIN%**" in sheet.row_values(row_index, 0, 1):
            start_data = row_index + 2
            break
    return start_data



def generate_coordinates_babel(smiles):
    try:
        pybelmol = pybel.readstring('smi', smiles)
        pybelmol.make3D()
        sdf_string = pybelmol.write("sdf")
        return sdf_string
    except:
        return ''

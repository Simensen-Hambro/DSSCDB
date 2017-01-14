from django.shortcuts import render
from .models import Article, Molecule, Spectrum, Performance
from django.utils.timezone import datetime
import requests
import bibtexparser
from django.core.exceptions import ObjectDoesNotExist
from django.core.exceptions import ValidationError

from django.contrib.auth.models import User
import re
from decimal import Decimal, InvalidOperation
from xlrd import open_workbook

def get_DOI_metadata(doi, user):
    url = 'http://dx.doi.org/' + doi
    headers = {'accept': 'application/x-bibtex', 'style': 'bibtex'}
    response_bytes = requests.get(url, headers=headers).content
    response_string = response_bytes.decode("utf-8")
    bibtex_object = bibtexparser.loads(response_string)

    try:
        article = bibtex_object.entries[0]
    except IndexError:
        print(doi)
        print(bibtex_object)

    new_article = {
        'author': article.get('author'),
        'title': article.get('title'),
        'journal': article.get('journal'),
        'volume': article.get('volume'),
        'doi': article.get('doi'),
        'pages': article.get('pages'),
        'electronic_id': article.get('ID'),
        'issue_nr': article.get('number'),
        'keywords': article.get('keywords'),
        # MISSING: KEYWORDS
        'year': datetime(year=int(article.get('year')), month=1, day=1),
        'user': user}

    return Article.objects.create(**new_article)

def to_decimal(string):
    if not isinstance(string, float):
        rep = re.compile('(\-?\d*\.?\d+)')
        result = rep.search(string)
        if result:
            return Decimal(result.group(0))
        else:
            return None
    else:
        return string


def load_raw_data():
    from xlrd import open_workbook
    user = User.objects.all()
    user = user[0]
    # user = request.user
    # user = 1

    book = open_workbook('odd.xls')
    sheet = book.sheet_by_index(0)

    for row_index in range(sheet.nrows):
        if "**%BEGIN%**" in sheet.row_values(row_index, 0, 1):
            start_data = row_index +2
            break

    for row_index in range(start_data, sheet.nrows):
        row = sheet.row_values(row_index, 0, 24)
        article_doi = row[0]

        try:
            article = Article.objects.get(doi=article_doi)
        except ObjectDoesNotExist:
            article = get_DOI_metadata(article_doi, user)

        try:
            molecule_data = {'user': user, 'smiles': row[15], 'inchi': row[16], 'keywords': row[20]}

            try:
                molecule = Molecule.objects.get(smiles=molecule_data['smiles'])
            except ObjectDoesNotExist:
                molecule = Molecule.objects.create(**molecule_data)

            spectum_data = {
                'absorption_maxima': to_decimal(row[17]), 'emission_maxima': to_decimal(row[18]), 'solvent': row[19],
                'molecule': molecule, 'user': user, 'article': article,
            }

            try:
                spectrum = Spectrum.objects.get(molecule=spectum_data['molecule'], user=user, article=article)
            except ObjectDoesNotExist:
                spectrum = Spectrum.objects.create(**spectum_data)

            performance_data = {
                'voc': to_decimal(row[1]),
                'jsc': to_decimal(row[2]),
                'ff': to_decimal(row[3]),
                'pce': to_decimal(row[4]),
                'electrolyte': row[5],
                'active_area': row[6],
                'co_adsorbent': row[7],
                'co_sensitizer': row[8],
                'semiconductor': row[9],
                'dye_loading': row[10],
                'exposure_time': row[11],
                'solar_simulator': row[12],
                'keywords': row[13],
                'comment': row[14],
                'molecule': molecule, 'user': user, 'article': article
            }

            Performance.objects.create(**performance_data)

        except IndexError:
            pass
            print("Error at line {}".format(row_index))
            # Incorrect formatting. Invisible linebreak.
            # Missing column
        except ValidationError:
            print("Error at line {}".format(row_index))
            # Incorrect formatting, illegal data format
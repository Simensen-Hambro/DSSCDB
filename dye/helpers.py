import re

import bibtexparser
import requests
from django.core.exceptions import ObjectDoesNotExist
from django.utils.timezone import datetime

from .models import Article


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
    }

    return Article.objects.create(**new_article)


def to_decimal(string):
    if not isinstance(string, float):
        illegal_characters = re.search('([^0-9^.^,^-])', string)
        if illegal_characters:
            return illegal_characters
        else:
            rep = re.compile('(\-?\d*\.?\d+)')

            result = rep.search(string)
            if result:
                return result.group(0)
            else:
                return None
    else:
        return string


def locate_start_data(sheet):
    # Locate control tag start row
    start_data = None
    for row_index in range(sheet.nrows):
        if "**%BEGIN%**" in sheet.row_values(row_index, 0, 1):
            start_data = row_index + 2
            break
    return start_data


def get_or_create_article(article_doi):
    try:
        article = Article.objects.get(doi__iexact=article_doi)
    except ObjectDoesNotExist:
        article = get_DOI_metadata(article_doi)
    return article


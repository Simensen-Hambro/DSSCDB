from django.shortcuts import render
from .models import Article, Molecule, Spectrum, Performance
from django.utils.timezone import datetime
import requests
import bibtexparser
from django.core.exceptions import ObjectDoesNotExist
from django.core.exceptions import ValidationError
from django.core.exceptions import FieldError

from django.contrib import messages
from django.shortcuts import reverse, render
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect

from .forms import ArticleForm, MoleculeForm, SpectrumForm, PerformanceForm
from django.contrib.auth.models import User
import re
from decimal import Decimal, InvalidOperation


def get_DOI_metadata(doi, user):
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


def locate_start_data(sheet):
    # Locate control tag start row
    start_data = -1
    for row_index in range(sheet.nrows):
        if "**%BEGIN%**" in sheet.row_values(row_index, 0, 1):
            start_data = row_index + 2
            break
    return start_data


def get_or_create_article(article_doi, user):
    try:
        article = Article.objects.get(doi=article_doi)
    except ObjectDoesNotExist:
        article = get_DOI_metadata(article_doi, user)
    return article


def get_or_create_molecule(molecule_data):
    try:
        molecule = Molecule.objects.get(smiles=molecule_data['smiles'])
    except ObjectDoesNotExist:
        molecule = Molecule.objects.create(**molecule_data)
    return molecule


def load_raw_data():
    from xlrd import open_workbook
    user = User.objects.all()
    user = user[0]
    # user = request.user
    # user = 1

    book = open_workbook('odd.xls')
    sheet = book.sheet_by_index(0)

    start_data = locate_start_data(sheet)

    if start_data == -1:
        messages.add_message(messages.ERROR, 'Could not find start-tag. Compare your sheet with the sample sheet.')

    for row_index in range(start_data, sheet.nrows):
        row = sheet.row_values(row_index, 0, 24)

        article_doi = row[0]
        article = get_or_create_article(article_doi, user)

        try:
            molecule_data = {'user': user, 'smiles': row[15], 'inchi': row[16], 'keywords': row[20]}
            molecule = get_or_create_molecule(**molecule_data)
            spectum_data = {
                'absorption_maxima': to_decimal(row[17]), 'emission_maxima': to_decimal(row[18]), 'solvent': row[19],
                'molecule': molecule, 'user': user, 'article': article,
            }

            try:
                spectrum = Spectrum.objects.get(molecule=spectum_data['molecule'], article=article)
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


@login_required
def upload_data(request):
    article_form = ArticleForm(request.POST or None)
    molecule_form = MoleculeForm(request.POST or None)
    spectrum_form = SpectrumForm(request.POST or None)
    performance_form = PerformanceForm(request.POST or None)

    article, molecule, spectrum, performance = None, None, None, None

    if request.method == "POST":
        if article_form.is_valid():
            try:
                article = get_or_create_article(article_form.cleaned_data.get('doi'), request.user)
                if not article:
                    article_form.add_error('doi', 'DOI not found')
                else:
                    try:
                        molecule_form.is_valid()
                        molecule = Molecule.objects.get(smiles=molecule_form.cleaned_data.get('smiles'))
                    except ObjectDoesNotExist:
                        if not molecule_form.is_valid():
                            raise FieldError

                    molecule = molecule_form.save(commit=False)
                    molecule.article = article
                    molecule.user = request.user

                    try:
                        spectrum = Spectrum.objects.get(molecule=molecule, article=article)
                    except ObjectDoesNotExist:
                        if not spectrum_form.is_valid():
                            raise FieldError

                    spectrum = spectrum_form.save(commit=False)
                    spectrum.molecule, spectrum.user, spectrum.article = molecule, request.user, article

                    try:
                        performance = Performance.objects.get(user=request.user, article=article, molecule=molecule)
                    except ObjectDoesNotExist:
                        if not performance_form.is_valid():
                            raise FieldError

                    performance = performance_form.save(commit=False)
                    performance.user, performance.article, performance.molecule = request.user, article, molecule

                    article.save()
                    molecule.save()
                    spectrum.save()
                    performance.save()
                    return redirect(reverse("dye:upload"))
            except FieldError:
                pass

    context = {
        'article': article,
        'molecule': molecule,
        'spectrum': spectrum,
        'performance': performance,
        'molecule_form': molecule_form,
        'article_form': article_form,
        'spectrum_form': spectrum_form,
        'performance_form': performance_form,
    }

    return render(request, 'dye/upload.html', context)

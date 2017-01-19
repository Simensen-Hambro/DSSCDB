from .models import Article, Molecule, Spectrum, Performance, Contribution
from django.utils.timezone import datetime
import requests
import bibtexparser
from django.core.exceptions import ObjectDoesNotExist
from django.core.exceptions import FieldError

from django.contrib import messages
from django.shortcuts import reverse, render
from django.contrib.auth.decorators import login_required
from django.shortcuts import redirect

from .forms import ArticleForm, MoleculeForm, SpectrumForm, PerformanceForm, SpreadsheetForm
import re
from django.db import transaction


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


def validate_raw_data(article_form, molecule_form, spectrum_form, performance_form, user):
    try:
        article = get_or_create_article(article_form.cleaned_data.get('doi'))
        if not article:
            article_form.add_error('doi', 'DOI not found')
            raise FieldError
        else:
            # article exists or was created
            try:
                # We try to fetch the molecule
                molecule_form.is_valid()
                molecule = Molecule.objects.get(smiles=molecule_form.data.get('smiles'))
            except ObjectDoesNotExist:
                # Molecule was not found, and therefore is_valid should now pass unless the user has erred
                if not molecule_form.is_valid():
                    raise FieldError
                molecule = molecule_form.save(commit=False)

            try:
                # Try to get spectrum
                spectrum = Spectrum.objects.get(molecule=molecule, article=article)
            except ObjectDoesNotExist:
                # Does not exist
                if not spectrum_form.is_valid():
                    raise FieldError
                spectrum = spectrum_form.save(commit=False)

            try:
                # Try to get performance
                # TODO: Does this "article, molecule" restriction make sense?
                performance = Performance.objects.get(article=article, molecule=molecule)
            except ObjectDoesNotExist:
                if not performance_form.is_valid():
                    raise FieldError

                performance = performance_form.save(commit=False)
                performance.user = user

            return {'article': article, 'molecule': molecule, 'spectrum': spectrum,
                    'performance': performance}

    except FieldError:
        return {'article': article_form, 'molecule': molecule_form, 'spectra': spectrum_form,
                'performance': performance_form}


@login_required
def single_upload(request):
    article_form = ArticleForm(request.POST or None)
    molecule_form = MoleculeForm(request.POST or None)
    spectrum_form = SpectrumForm(request.POST or None)
    performance_form = PerformanceForm(request.POST or None)
    forms = {'article_form': article_form, 'molecule_form': molecule_form, 'spectrum_form': spectrum_form,
             'performance_form': performance_form}

    if request.method == "POST":
        if article_form.is_valid():
            data_objects = validate_raw_data(user=request.user, **forms)
            Contribution.objects.create_from_data([data_objects], user=request.user)
            if isinstance(data_objects.get('article'), Article):
                messages.add_message(request, messages.SUCCESS,
                                     'The data was uploaded and is awaiting review. Thank you!')
                return redirect(reverse("dye:single-upload"))

    context = {
        'molecule_form': molecule_form,
        'article_form': article_form,
        'spectrum_form': spectrum_form,
        'performance_form': performance_form,
    }

    return render(request, 'dye/single-upload.html', context)


@login_required
def file_upload(request):
    user = request.user
    file_form = SpreadsheetForm(request.POST or None, request.FILES or None)
    if request.method == 'POST':
        if file_form.is_valid():
            # Posted valid data
            from xlrd import open_workbook
            upload = file_form.save(commit=False)
            upload.user = user
            upload.save()

            book = open_workbook('./' + str(upload.file))
            sheet = book.sheet_by_index(0)

            start_data = locate_start_data(sheet)
            if not start_data:
                messages.add_message(messages.ERROR,
                                     'Could not find start-tag. Compare your sheet with the sample sheet.')
                return redirect(reverse('dye:file-upload'))

            results = []
            for row_index in range(start_data, sheet.nrows):
                row = sheet.row_values(row_index, 0, 24)

                try:
                    # Populate article, molecule, spectrum and performance forms with the data from the user

                    article_form = ArticleForm({'doi': row[0]})
                    article_form.is_valid()
                    molecule_form = MoleculeForm(
                        {'user': user, 'smiles': row[15], 'inchi': row[16], 'keywords': row[20]})
                    spectrum_form = SpectrumForm({
                        'absorption_maxima': to_decimal(row[17]), 'emission_maxima': to_decimal(row[18]),
                        'solvent': row[19]
                    })
                    performance_form = PerformanceForm({
                        'voc': to_decimal(row[1]), 'jsc': to_decimal(row[2]), 'ff': to_decimal(row[3]),
                        'pce': to_decimal(row[4]), 'electrolyte': row[5], 'active_area': row[6], 'co_adsorbent': row[7],
                        'co_sensitizer': row[8], 'semiconductor': row[9], 'dye_loading': row[10],
                        'exposure_time': row[11], 'solar_simulator': row[12], 'keywords': row[13], 'comment': row[14]
                    })

                    forms = {'article_form': article_form, 'molecule_form': molecule_form,
                             'spectrum_form': spectrum_form,
                             'performance_form': performance_form}

                    results.append(validate_raw_data(user=user, **forms))

                except IndexError:
                    # Failed to get some value, raise error.
                    messages.add_message(request, messages.ERROR,
                                         'Critical error at row {}. '.format(row_index))
            # Iterate over attribute error, for every row
            errors = []
            for row_nr, row_data in enumerate(results):
                for instance in row_data:
                    if not hasattr(instance, 'pk'):
                        for k, v in instance.errors.items():
                            errors.append(
                                {'row': start_data + row_nr, 'key': k.replace('_', ' ').title(), 'message': v})

            if errors:
                messages.add_message(request, messages.ERROR, 'Upload failed')
                return render(request, 'dye/file-upload.html', context={'file_form': file_form, 'errors': errors})
            else:
                Contribution.objects.create(results)
                messages.add_message(request, messages.SUCCESS,
                                     'The data was uploaded and is awaiting review. Thank you!')
            return redirect(reverse('dye:file-upload'))

    return render(request, 'dye/file-upload.html', context={'file_form': file_form})

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import FieldError
from django.core.exceptions import ObjectDoesNotExist
from django.shortcuts import redirect
from django.shortcuts import reverse, render
from django.forms import modelformset_factory
from django.shortcuts import reverse, render, Http404
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger


from .forms import ArticleForm, MoleculeForm, SpectrumForm, PerformanceForm, SpreadsheetForm, ApprovalForm
from .helpers import get_or_create_article, locate_start_data, to_decimal
from .models import Article, Molecule, Spectrum, Performance, Contribution, APPROVAL_STATES

from django.db import transaction


@transaction.atomic
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
                molecule = Molecule.objects.get(inchi=molecule_form.data.get('inchi'))
            except ObjectDoesNotExist:
                # Molecule was not found, and therefore is_valid should now pass unless the user has erred
                if not molecule_form.is_valid():
                    raise FieldError
                molecule_form.article = article
                molecule = molecule_form.save()

            try:
                # Try to get spectrum
                spectrum = Spectrum.objects.get(molecule=molecule, article=article)
            except ObjectDoesNotExist:
                # Does not exist
                if not spectrum_form.is_valid():
                    raise FieldError
                spectrum = spectrum_form.save(commit=False)
                spectrum.article, spectrum.molecule = article, molecule
                spectrum.save()
            try:
                # Try to get performance
                # TODO: Does this "article, molecule" restriction make sense?
                performance = Performance.objects.get(article=article, molecule=molecule)
            except ObjectDoesNotExist:
                if not performance_form.is_valid():
                    raise FieldError

                performance = performance_form.save(commit=False)
                performance.article, performance.molecule, performance.user = article, molecule, user
                performance.save()

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

            book = open_workbook("media/" + str(upload.file))
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
                for instance in row_data.values():
                    if not hasattr(instance, 'pk'):
                        for k, v in instance.errors.items():
                            errors.append(
                                {'row': start_data + 1 + row_nr, 'key': k.replace('_', ' ').title(), 'message': v})

            if errors:
                messages.add_message(request, messages.ERROR, 'Upload failed')
                return render(request, 'dye/file-upload.html', context={'file_form': file_form, 'errors': errors})
            else:
                Contribution.objects.create_from_data(results, user=user)
                messages.add_message(request, messages.SUCCESS,
                                     'The data was uploaded and is awaiting review. Thank you!')
            return redirect(reverse('dye:file-upload'))

    return render(request, 'dye/file-upload.html', context={'file_form': file_form})

@login_required
def performance_list(request):
    from itertools import chain
    performance_list = Performance.objects.all()
    for i in range(10):
        performance_list = list(chain(performance_list, Performance.objects.all()))
    paginator = Paginator(performance_list, 10)

    page = request.GET.get('page')
    try:
        performances = paginator.page(page)
        page = int(page)
    except PageNotAnInteger:
        performances = paginator.page(1)
        page = 1
    except EmptyPage:
        performances = paginator.page(paginator.num_pages)
        page = paginator.num_pages

    if page - 5 < 1:
        first_page = 1
    else:
        first_page = page - 5

    if page + 5 > paginator.num_pages:
        last_page = paginator.num_pages
    else:
        last_page = page + 5

    context = {
        'performances': performances,
        'pages': [i for i in range(first_page, last_page + 1)],
        'num_pages': paginator.num_pages,
    }

    return render(request, 'dye/performance_list.html', context)

@login_required
def performance_details(request, short_id):
    try:
        performance = Performance.objects.get(short_id=short_id)
    except Performance.DoesNotExist:
        raise Http404

    context = {
        'performance': performance
    }

    return render(request, 'dye/performance_detail.html', context)

@login_required
def contributions_evaluation_overview(request):
    to_evaluate = Contribution.objects.filter(status__in=[APPROVAL_STATES.DENIED, APPROVAL_STATES.WAITING])
    ApprovalFormSet = modelformset_factory(Contribution, fields=('status', 'molecules'))

    if request.method == 'POST':
        formset = ApprovalFormSet(request.POST)
        if formset.is_valid():
            instances = formset.save(commit=False)
            for instance in instances.changed_objects():
                instance.save()

    else:
        formset = ApprovalFormSet(
            queryset=Contribution.objects.filter(status__in=[APPROVAL_STATES.DENIED, APPROVAL_STATES.WAITING]))
    context = {
        'contributions': to_evaluate,
        'formset': formset,
    }
    return render(request, 'dye/evaluate_contributions.html', context=context)

@login_required
def single_contribution_evaluation(request, contribution):
    contribution = Contribution.objects.get(pk=contribution)
    performances = contribution.performances.all()

    approval_form = ApprovalForm(request.POST or None, instance=contribution)

    if request.method == 'POST':
        if approval_form.is_valid():
            approval_form.save()
            messages.add_message(request, messages.SUCCESS,
                                 'The contribution has been marked as {}'.format(
                                     APPROVAL_STATES.for_value(contribution.status).display))
            return redirect(reverse("dye:evaluate-contributions"))
    return render(request, 'dye/single_evaluation.html', context={'approval_form': approval_form})

@login_required
def my_contributions(request):
    my_contributions = Contribution.objects.filter(user=request.user)
    return render(request, 'dye/my_contributions.html', context={'contributions':my_contributions})


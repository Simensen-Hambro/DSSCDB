from django.conf import settings
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import FieldError
from django.core.exceptions import ObjectDoesNotExist
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.db import IntegrityError, transaction
from django.forms import modelformset_factory
from django.shortcuts import reverse, render, redirect, Http404

from .forms import *
from .helpers import get_or_create_article, locate_start_data, to_decimal
from .models import Molecule, Spectrum, Performance, Contribution, APPROVAL_STATES


def validate_raw_data(article_form, molecule_form, spectrum_form, performance_form, user):
    """
    :param *form: Forms necessary to create a "performance" object  
    :return: Tuple: (Passed, [(object_1, created_1), (object_2, created_2), ... ]).
        Passed: True if there were no errors during the process, false otherwise
        Objects: list of tuples containing the object it self as well as a boolean if it was created now
    """
    try:
        article, created_article = get_or_create_article(article_form.cleaned_data.get('doi'))
        if not article:
            article_form.add_error('doi', 'DOI not found')
            raise FieldError
        else:
            # article exists or was created
            created_molecule, created_spectrum, created_performance = False, False, False
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
                created_molecule = True

            try:
                # TODO: Remove one to one constraint on spectrum.
                # Try to get spectrum
                spectrum = Spectrum.objects.get(molecule=molecule)
            except ObjectDoesNotExist:
                # Does not exist
                if not spectrum_form.is_valid():
                    raise FieldError
                spectrum = spectrum_form.save(commit=False)
                spectrum.article, spectrum.molecule = article, molecule
                spectrum.save()
                created_spectrum = True
            try:
                # Try to get performance
                performance = Performance.objects.get(article=article, molecule=molecule,
                                                      voc=performance_form.data.get('voc'),
                                                      jsc=performance_form.data.get('jsc'),
                                                      ff=performance_form.data.get('ff'),
                                                      pce=performance_form.data.get('pce'))
            except ObjectDoesNotExist:
                if not performance_form.is_valid():
                    raise FieldError

                performance = performance_form.save(commit=False)
                performance.article, performance.molecule, performance.user = article, molecule, user
                performance.save()
                created_performance = True

            return True, [(article, created_article), (molecule, created_molecule), (spectrum, created_spectrum),
                          (performance, created_performance)]

    except FieldError:
        return False, [(article_form, False), (molecule_form, False), (spectrum_form, False), (performance_form, False)]


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
            passed, data_objects = validate_raw_data(user=request.user, **forms)
            Contribution.objects.create_from_data([data_objects], user=request.user)

            # If FieldError in validate_raw_data, the article is an article_form
            if passed is True:
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
            from xlrd import open_workbook, XLRDError

            upload = file_form.save(commit=False)
            upload.user = user
            upload.save()

            try:
                book = open_workbook(settings.MEDIA_ROOT + '/' + str(upload.file))
            except XLRDError:
                messages.add_message(request, messages.ERROR,
                                     'The file was not recognized as a valid spreadsheet file. '
                                     'Please download the sample file and try again.')
                return redirect(reverse('dye:file-upload'))

            sheet = book.sheet_by_index(0)

            start_data = locate_start_data(sheet)
            if not start_data:
                messages.add_message(request, messages.ERROR,
                                     'Could not find start-tag. Compare your sheet with the sample sheet.')
                return redirect(reverse('dye:file-upload'))

            results, total_status = [], True

            try:
                with transaction.atomic():
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
                                'pce': to_decimal(row[4]), 'electrolyte': row[5], 'active_area': row[6],
                                'co_adsorbent': row[7],
                                'co_sensitizer': row[8], 'semiconductor': row[9], 'dye_loading': row[10],
                                'exposure_time': row[11], 'solar_simulator': row[12], 'keywords': row[13],
                                'comment': row[14]
                            })

                            forms = {'article_form': article_form, 'molecule_form': molecule_form,
                                     'spectrum_form': spectrum_form,
                                     'performance_form': performance_form}

                            passed, data = validate_raw_data(user=user, **forms)
                            _, _, _, r = data
                            _, c = r
                            if not c:
                                print("Row {} did not yield performance".format(row_index))
                            results.append(data)

                            if passed is not True:
                                total_status = False

                        except IndexError:
                            # Failed to get some value, raise error.
                            total_status = False
                            messages.add_message(request, messages.ERROR,
                                                 'Critical error at row {}. '.format(row_index))

                    if total_status is not True:
                        raise IntegrityError
            except IntegrityError:
                pass

            if total_status is not True:
                # Iterate over attribute error, for every row
                errors = []
                # _, objects = zip(*results)
                for row_nr, row_data in enumerate(results):
                    row_data, _ = zip(*row_data)
                    for instance in row_data:
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
def contributions_evaluation_overview(request):
    to_evaluate = Contribution.objects.filter(status__in=[APPROVAL_STATES.DENIED, APPROVAL_STATES.WAITING])
    ApprovalFormSet = modelformset_factory(Contribution, fields=('status',))

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
        'admin': True,
    }
    return render(request, 'dye/evaluate_contributions.html', context=context)


@login_required
def contribution_performances(request, short_id):
    contribution = Contribution.objects.get(short_id=short_id)
    atomic_contrib_performances = contribution.items.filter(content_type=ContentType.objects.get_for_model(Performance))
    performances = Performance.objects.filter(contribution__in=atomic_contrib_performances)
    approval_form = None

    if request.user.has_perm('dye.contribution.set_contribution_status'):
        approval_form = ApprovalForm(request.POST or None, instance=contribution)

        if request.method == 'POST':
            if approval_form.is_valid():
                approval_form.save()
                messages.add_message(request, messages.SUCCESS,
                                     'The contribution has been marked as {}'.format(
                                         APPROVAL_STATES.for_value(contribution.status).display))
                return redirect(reverse("dye:evaluate-contributions"))

    return render(request, 'dye/single_evaluation.html',
                  context={'approval_form': approval_form, 'performances': performances})


@login_required
def my_contributions(request):
    contributions = Contribution.objects.filter(user=request.user)
    return render(request, 'dye/my_contributions.html', context={'contributions': contributions})


@login_required
def performance_list(request):
    context = paginate_performances(request, get_performances(), {})
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
def performance_search(request):
    context = {}
    context['search'] = False
    if request.method == 'POST':
        range_form = PerformanceRangeSearchForm(request.POST)
        keyword_form = PerformanceKeywordSearchForm(request.POST)
        structure_form = PerformanceStructureSearchForm(request.POST)
        if range_form.is_valid() and keyword_form.is_valid() and structure_form.is_valid():
            performances = get_performances(**range_form.cleaned_data, **keyword_form.cleaned_data,
                                            **structure_form.cleaned_data)
            context = paginate_performances(request, performances, context)
            context['hits'] = performances.count()
            context['search'] = True
    else:
        range_form = PerformanceRangeSearchForm()
        keyword_form = PerformanceKeywordSearchForm()
        structure_form = PerformanceStructureSearchForm()
        performances = get_performances()
        context = paginate_performances(request, performances, context)

    context['range_form'] = range_form
    context['keyword_form'] = keyword_form
    context['structure_form'] = structure_form

    return render(request, 'dye/performance_list.html', context)


def get_performances(**search):
    if search.get('status'):
        # performances = Performance.objects.filter(status=search.get('status'))
        pass
    else:
        pass
        # performances = Performance.objects.filter(status=APPROVAL_STATES.APPROVED)

    performances = Performance.objects.filter(status=APPROVAL_STATES.APPROVED)

    # Search after keyword
    if search.get('keyword'):
        performances = performances.filter(keywords__icontains=search.get('keyword'))

    # Structure search
    if search.get('smiles'):
        performances = performances.filter(molecule__smiles__icontains=search.get('smiles'))

    # Search after different range criterias
    if search.get('min_voc'):
        performances = performances.filter(voc__gte=search.get('min_voc'))
    if search.get('max_voc'):
        performances = performances.filter(voc__lte=search.get('max_voc'))
    if search.get('min_jsc'):
        performances = performances.filter(jsc__gte=search.get('min_jsc'))
    if search.get('max_jsc'):
        performances = performances.filter(jsc__lte=search.get('max_jsc'))
    if search.get('min_ff'):
        performances = performances.filter(ff__gte=search.get('min_ff'))
    if search.get('max_ff'):
        performances = performances.filter(ff__lte=search.get('max_ff'))
    if search.get('min_pce'):
        performances = performances.filter(pce__lte=search.get('min_pce'))
    if search.get('max_pce'):
        performances = performances.filter(pce__lte=search.get('max_pce'))

    return performances


def paginate_performances(request, performance_list, context):
    """
    For a given set of performances returns the context with pagination
    """
    paginator = Paginator(performance_list, settings.PAGINATION_NUMBER)

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

    context['performances'] = performances
    context['pages'] = [i for i in range(first_page, last_page + 1)]
    context['num_pages'] = paginator.num_pages

    return context

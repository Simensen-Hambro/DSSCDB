from django import forms

from .helpers import get_DOI_metadata
from .models import Molecule, Spectrum, Performance, Spreadsheet, Contribution, Article


class ArticleForm(forms.Form):
    doi = forms.CharField(max_length=500, label="DOI", required=True)

    def is_valid(self):
        super().is_valid()

        article_doi = self.cleaned_data.get('doi')
        article = Article.objects.filter(doi__iexact=article_doi).first()

        if article:
            self.model_instance = article
            self.created = False
            return True
        else:
            # TODO: Add try and expect if connection to DOI is not found
            article_data = get_DOI_metadata(article_doi)
            if not article_data:
                self.add_error('doi', 'DOI not found')
                return False

            article_model = ArticleModelForm(article_data)
            if article_model.is_valid():
                article_model.save()
                self.model_instance = article_model.instance
                self.created = True
                return True
            else:
                erred_fields = ', '.join(article_model.errors.keys())
                self.add_error('doi',
                               'DOI provided has incomplete ({}) data. Please contact us regarding this.'.format(
                                   erred_fields))
                return False
                # except TypeError:
                #    self.add_error('doi', 'DOI not found')

    def get_model(self):
        return self.model_instance, self.created


class ArticleModelForm(forms.ModelForm):
    class Meta:
        model = Article
        fields = [
            'author',
            'title',
            'journal',
            'volume',
            'doi',
            'pages',
            'issue_nr',
            'eid',
            'year',
            'electronic_id',
            'keywords',
        ]


class MoleculeForm(forms.ModelForm):
    class Meta:
        model = Molecule
        fields = [
            'smiles',
            'inchi',
            'keywords'
        ]

    def validate_unique(self):
        pass


class SpectrumForm(forms.ModelForm):
    class Meta:
        model = Spectrum
        fields = [
            'absorption_maxima',
            'emission_maxima',
            'solvent',
        ]


class PerformanceForm(forms.ModelForm):
    class Meta:
        model = Performance
        fields = [
            'voc',
            'jsc',
            'ff',
            'pce',
            'electrolyte',
            'active_area',
            'co_adsorbent',
            'co_sensitizer',
            'semiconductor',
            'dye_loading',
            'exposure_time',
            'solar_simulator',
            'comment',
        ]

    def is_valid(self, article, molecule):
        super().is_valid()
        if self.errors:
            return False

        # Try to get existing performance, add error if duplicate found
        try:
            performance = Performance.objects.get(article=article, molecule=molecule,
                                                  voc=str(self.data.get('voc')),
                                                  jsc=str(self.data.get('jsc')),
                                                  ff=str(self.data.get('ff')),
                                                  pce=str(self.data.get('pce')))
            if performance:
                self.add_error('voc', 'The performance measure exists already for the given molecule')
                return False
        except Performance.DoesNotExist:
            # All ok
            pass

        return True

    def clean(self):
        super().clean()


class SpreadsheetForm(forms.ModelForm):
    class Meta:
        model = Spreadsheet
        fields = [
            'file'
        ]


class ApprovalForm(forms.ModelForm):
    confirm = forms.BooleanField(required=True, label="I've thoroughly checked the data", )

    class Meta:
        model = Contribution
        fields = [
            'status'
        ]


class PerformanceRangeSearchForm(forms.Form):
    min_voc = forms.DecimalField(label='min VOC', decimal_places=4, max_digits=15, required=False)
    max_voc = forms.DecimalField(label='max VOC', decimal_places=4, max_digits=15, required=False)
    min_jsc = forms.DecimalField(label='min JSC', decimal_places=5, max_digits=15, required=False)
    max_jsc = forms.DecimalField(label='max JSC', decimal_places=5, max_digits=15, required=False)
    min_ff = forms.DecimalField(label='min FF', decimal_places=5, max_digits=13, required=False)
    max_ff = forms.DecimalField(label='max FF', decimal_places=5, max_digits=13, required=False)
    min_pce = forms.DecimalField(label='min PCE', decimal_places=5, max_digits=13, required=False)
    max_pce = forms.DecimalField(label='max PCE', decimal_places=5, max_digits=13, required=False)


class PerformanceKeywordSearchForm(forms.Form):
    keyword = forms.CharField(max_length=1000,
                              required=False,
                              widget=forms.TextInput(attrs={'placeholder': 'Free text search'}))


class PerformanceStructureSearchForm(forms.Form):
    smiles = forms.CharField(max_length=1000, required=False)
    complete_molecule = forms.BooleanField(required=False)

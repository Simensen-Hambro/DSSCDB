from django import forms

from .models import Molecule, Spectrum, Performance, Spreadsheet, Contribution


class ArticleForm(forms.Form):
    doi = forms.CharField(max_length=500, label="DOI")


class MoleculeForm(forms.ModelForm):
    class Meta:
        model = Molecule
        fields = [
            'smiles',
            'inchi',
            'keywords'
        ]


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
            'keywords',
            'comment',
        ]


class SpreadsheetForm(forms.ModelForm):
    class Meta:
        model = Spreadsheet
        fields = [
            'file'
        ]


class ApprovalForm(forms.ModelForm):
    confirm = forms.BooleanField(required=True, label="I've thoroughly checked the data")

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

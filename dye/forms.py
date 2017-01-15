from django import forms
from .models import Molecule, Spectrum, Performance, Spreadsheet


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

from django.db import models
from sorl.thumbnail import ImageField
from django.contrib.auth.models import User


class Molecule(models.Model):
    smiles = models.CharField(max_length=1000, verbose_name='SMILES')
    inchi = models.CharField(max_length=1000, verbose_name='INCHI')
    image = ImageField(upload_to='molecules', verbose_name='Picture')
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)
    user = models.ForeignKey(User, related_name='molecules')


class Article(models.Model):
    author = models.CharField(max_length=1000)
    title = models.CharField(max_length=1000)
    journal = models.CharField(max_length=250)
    volume = models.PositiveIntegerField()
    doi = models.CharField(max_length=500, verbose_name='DOI')
    pages = models.CharField(max_length=20)
    issue_nr = models.PositiveIntegerField()
    eid = models.PositiveIntegerField(blank=True, null=True)
    year = models.DateField()
    electronic_id = models.CharField(max_length=250)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)
    user = models.ForeignKey(User, related_name='articles')


class Spectrum(models.Model):
    absorption_maxima = models.FloatField()
    emission_maxima = models.FloatField()
    solvent = models.CharField(max_length=100)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    molecule = models.ForeignKey(Molecule)
    user = models.ForeignKey(User, related_name='spectra')
    article = models.ForeignKey(Article, related_name='spectra')

    class Meta:
        verbose_name = "Molecule spectrum"
        verbose_name_plural = "Molecule spectra"


class Performance(models.Model):
    pce = models.FloatField(verbose_name='PCE')
    jsc = models.FloatField(verbose_name='JSC')
    voc = models.FloatField(verbose_name='VOC')
    ff = models.FloatField(verbose_name='FF')
    electrolyte = models.CharField(max_length=1000)
    active_area = models.CharField(max_length=30)
    co_adsorbent = models.CharField(max_length=250)
    co_sensitizer = models.CharField(max_length=1000)
    semiconductor = models.CharField(max_length=1000)
    dye_loading = models.CharField(max_length=1000)
    exposure_time = models.CharField(max_length=500)
    solar_simulator = models.CharField(max_length=1000)

    comment = models.TextField()
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)
    user = models.ForeignKey(User, related_name='performances')
    article = models.ForeignKey(Article, related_name='performances')
    molecule = models.ForeignKey(Molecule, related_name='performances')

    class Meta:
        verbose_name = "DSSC performance"

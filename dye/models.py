from django.db import models
from sorl.thumbnail import ImageField
from django.contrib.auth.models import User


class Molecule(models.Model):
    smiles = models.CharField(max_length=1000, verbose_name='SMILES', blank=True, null=True)
    # See Record 177, 178 for this "blank"
    inchi = models.CharField(max_length=1000, verbose_name='INCHI', blank=True, null=True)
    image = ImageField(upload_to='molecules', verbose_name='Picture')
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)
    user = models.ForeignKey(User, related_name='molecules')

    def __str__(self):
        return self.inchi


class Article(models.Model):
    author = models.CharField(max_length=1000)
    title = models.CharField(max_length=1000)
    journal = models.CharField(max_length=250)
    volume = models.CharField(max_length=100)
    doi = models.CharField(max_length=500, verbose_name='DOI')
    pages = models.CharField(max_length=20)
    issue_nr = models.CharField(max_length=100, blank=True, null=True)
    eid = models.CharField(blank=True, null=True, max_length=100)
    year = models.DateField()
    electronic_id = models.CharField(max_length=250)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)
    user = models.ForeignKey(User, related_name='articles')

    def __str__(self):
        return '[{}] - "{}", vol. {}, issue {}, {} '.format(self.doi, self.title, self.volume, self.issue_nr, self.year)


class Spectrum(models.Model):
    absorption_maxima = models.DecimalField(blank=True, null=True, decimal_places=4, max_digits=10)
    emission_maxima = models.DecimalField(blank=True, null=True, decimal_places=4, max_digits=10)
    solvent = models.CharField(max_length=100)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    molecule = models.ForeignKey(Molecule)
    user = models.ForeignKey(User, related_name='spectra')
    article = models.ForeignKey(Article, related_name='spectra')

    class Meta:
        verbose_name = "Molecule spectrum"
        verbose_name_plural = "Molecule spectra"

    def __str__(self):
        return '{} - abs. max {}, emi. max {}'.format(self.molecule, self.absorption_maxima, self.emission_maxima)


class Performance(models.Model):
    voc = models.DecimalField(verbose_name='VOC', decimal_places=4, max_digits=15)
    jsc = models.DecimalField(verbose_name='JSC', decimal_places=10, max_digits=15)
    ff = models.DecimalField(verbose_name='FF', decimal_places=10, max_digits=13)
    pce = models.DecimalField(verbose_name='PCE', decimal_places=10, max_digits=13)
    electrolyte = models.CharField(max_length=1000)
    active_area = models.CharField(max_length=30)
    co_adsorbent = models.CharField(max_length=250, blank=True, null=True)
    co_sensitizer = models.CharField(max_length=1000, blank=True, null=True)
    semiconductor = models.CharField(max_length=1000)
    dye_loading = models.CharField(max_length=1000)
    exposure_time = models.CharField(max_length=500)
    solar_simulator = models.CharField(max_length=1000)
    keywords = models.CharField(max_length=1000, blank=True, null=True)

    comment = models.TextField()
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    user = models.ForeignKey(User, related_name='performances')
    article = models.ForeignKey(Article, related_name='performances')
    molecule = models.ForeignKey(Molecule, related_name='performances')

    def __str__(self):
        return str(self.molecule)

    class Meta:
        verbose_name = "DSSC performance"
        permissions = (
            ("upload_performance_data", "Can upload performance data"),
        )

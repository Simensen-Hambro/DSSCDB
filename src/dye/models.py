import uuid

from django.contrib.auth.models import User
from django.db import models
from django.shortcuts import reverse
from extended_choices import Choices
from sorl.thumbnail import ImageField
from tinyuuidfield.fields import TinyUUIDField

from .validators import validate_inchi, validate_smiles

APPROVAL_STATES = Choices(
    ('WAITING', 1, 'Waiting'),
    ('APPROVED', 2, 'Approved'),
    ('DENIED', 3, 'Denied'),
)

class Spreadsheet(models.Model):
    file = models.FileField(upload_to='spreadsheets')
    user = models.ForeignKey(User)
    created = models.DateTimeField(auto_now_add=True)


class Molecule(models.Model):
    smiles = models.CharField(max_length=1000, verbose_name='SMILES', unique=True, help_text="Example field help text.")
    inchi = models.CharField(max_length=1000, verbose_name='INCHI', unique=True)
    image = ImageField(upload_to='molecules', verbose_name='Picture', blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)

    def __str__(self):
        return self.inchi


class Article(models.Model):
    author = models.CharField(max_length=1000)
    title = models.CharField(max_length=1000)
    journal = models.CharField(max_length=250)
    volume = models.CharField(max_length=100)
    doi = models.CharField(max_length=500, verbose_name='DOI', unique=True)
    pages = models.CharField(max_length=20)
    issue_nr = models.CharField(max_length=100, blank=True, null=True)
    eid = models.CharField(blank=True, null=True, max_length=100)
    year = models.DateField()
    electronic_id = models.CharField(max_length=250)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)

    def __str__(self):
        return '{} - "{}", vol. {}, issue {}, {} '.format(self.doi, self.title, self.volume, self.issue_nr, self.year)


class Spectrum(models.Model):
    absorption_maxima = models.DecimalField(blank=True, null=True, decimal_places=4, max_digits=10, help_text="[kg/s]")
    emission_maxima = models.DecimalField(blank=True, null=True, decimal_places=4, max_digits=10)
    solvent = models.CharField(max_length=100)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    # Foreignkey, not one to one
    molecule = models.OneToOneField(Molecule, related_name='spectrum')
    article = models.ForeignKey(Article, related_name='spectra')

    status = models.PositiveSmallIntegerField(choices=APPROVAL_STATES, default=APPROVAL_STATES.WAITING)

    class Meta:
        unique_together = ('molecule', 'article')
        verbose_name = "Molecule spectrum"
        verbose_name_plural = "Molecule spectra"

    def set_status(self, status):
        self.status = status

    def __str__(self):
        return '{} - abs. max {}, emi. max {}'.format(self.molecule, self.absorption_maxima, self.emission_maxima)


class Performance(models.Model):
    voc = models.DecimalField(verbose_name='VOC', decimal_places=4, max_digits=15)
    jsc = models.DecimalField(verbose_name='JSC', decimal_places=5, max_digits=15)
    ff = models.DecimalField(verbose_name='FF', decimal_places=5, max_digits=13)
    pce = models.DecimalField(verbose_name='PCE', decimal_places=5, max_digits=13)
    electrolyte = models.CharField(max_length=1000)
    active_area = models.CharField(max_length=30, help_text='cm2')
    co_adsorbent = models.CharField(max_length=250, blank=True, null=True)
    co_sensitizer = models.CharField(max_length=1000, blank=True, null=True)
    semiconductor = models.CharField(max_length=1000)
    dye_loading = models.CharField(max_length=1000, help_text='nmol/cm2')
    exposure_time = models.CharField(max_length=500)
    solar_simulator = models.CharField(max_length=1000)
    keywords = models.CharField(max_length=1000, blank=True, null=True)

    short_id = TinyUUIDField(length=10)
    comment = models.TextField()
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    article = models.ForeignKey(Article, related_name='performances')
    molecule = models.ForeignKey(Molecule, related_name='performances')

    status = models.PositiveSmallIntegerField(choices=APPROVAL_STATES, default=APPROVAL_STATES.WAITING)

    def set_status(self, status):
        self.status = status

    def __str__(self):
        return str(self.molecule)

    def get_absolute_url(self):
        return reverse("dye:performance-detail", kwargs={'short_id': self.short_id})

    class Meta:
        verbose_name = "DSSC performance"


class ContributionManager(models.Manager):
    def create_from_data(self, data_entry, *args, **kwargs):
        contribution = self.create(*args, **kwargs)
        for row in data_entry:
            article, molecule, spectrum, performance = row.get('article'), row.get('molecule'), \
                                                       row.get('spectrum'), row.get('performance')
            contribution.articles.add(article)
            contribution.molecules.add(molecule)
            contribution.specta.add(spectrum)
            contribution.performances.add(performance)
        return contribution


class Contribution(models.Model):
    user = models.ForeignKey(User, related_name='contributions')
    performances = models.ManyToManyField(Performance, null=True, blank=True)
    articles = models.ManyToManyField(Article, null=True, blank=True)
    specta = models.ManyToManyField(Spectrum, null=True, blank=True)
    molecules = models.ManyToManyField(Molecule)

    short_id = TinyUUIDField(length=10)
    created = models.DateTimeField(auto_now_add=True)
    status = models.PositiveSmallIntegerField(choices=APPROVAL_STATES, default=APPROVAL_STATES.WAITING)

    objects = ContributionManager()

    class Meta:
        permissions = (
            ("upload_performance_data", "Can upload performance data"),
            ("set_contribution_status", "Can change status on contributions")
        )

    def approve(self):
        self.set_status(APPROVAL_STATES.APPROVED)

    def deny(self):
        self.set_status(APPROVAL_STATES.DENIED)

    def wait(self):
        self.set_status(APPROVAL_STATES.WAITING)

    def set_approval_state(self, status):
        for obj in [self.performances.all(), self.specta.all()]:
            for item in obj:
                item.set_status(status)
                item.save()

    def save(self, *args, **kwargs):
        super(Contribution, self).save(*args, **kwargs)
        if self.pk:
            self.set_approval_state(self.status)

    def get_absolute_url(self):
        return reverse("dye:single-evaluation", kwargs={'short_id': self.short_id})

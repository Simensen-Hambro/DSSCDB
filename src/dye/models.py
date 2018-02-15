from itertools import chain
from operator import itemgetter

import pybel
from django.contrib.auth.models import User
from django.contrib.contenttypes.fields import GenericForeignKey, GenericRelation
from django.contrib.contenttypes.models import ContentType
from django.db import models
from django.shortcuts import reverse
from extended_choices import Choices
from rdkit import Chem
from tinyuuidfield.fields import TinyUUIDField

from .validators import validate_smiles

APPROVAL_STATES = Choices(
    ('WAITING', 1, 'Waiting'),
    ('APPROVED', 2, 'Approved'),
    ('DENIED', 3, 'Denied'),
)


class AtomicContribution(models.Model):
    """
        A model that ties together any uploaded object type, with a Contribution object    
    """
    content_type = models.ForeignKey(ContentType, on_delete=models.CASCADE)
    object_id = models.PositiveIntegerField()
    content_object = GenericForeignKey('content_type', 'object_id')


class Data(models.Model):
    status = models.PositiveSmallIntegerField(choices=APPROVAL_STATES, default=APPROVAL_STATES.WAITING)
    upload = GenericRelation('AtomicContribution')

    class Meta:
        abstract = True

    def set_status(self, status):
        self.status = status


class Spreadsheet(Data):
    file = models.FileField(upload_to='spreadsheets')
    user = models.ForeignKey(User)
    created = models.DateTimeField(auto_now_add=True)


class MoleculeManager(models.Manager):
    def substructure_search(self, query_smiles):
        # Convert molecule SMILES by passing it through a RDKIT molecule object to ensure that it is kukulized
        query_molecule = pybel.readstring('smiles', query_smiles)
        query_canonical = query_molecule.write('can')
        query_smarts = pybel.Smarts(query_canonical)

        all_molecules = self.get_queryset().filter(status=APPROVAL_STATES.APPROVED)
        result = []
        for molecule in all_molecules:
            molecule_pybel = pybel.readstring('smiles', molecule.smiles)
            matches = query_smarts.findall(molecule_pybel)
            if matches:
                result.append(molecule.pk)

        return self.filter(pk__in=result)

    def similarity_search(self, query_smiles, tanimoto_threshold=0.8):
        # Fetch all molecules
        all_molecules = self.get_queryset().filter(status=APPROVAL_STATES.APPROVED)

        # The fingerprint bits are stored as comma separated string, split it, and int all values, for all molecules
        # Result is stored as a tuple with the primary key of the molecule
        all_bitlists = [([int(b) for b in m.fingerprint.split(',')], m.id) for m in all_molecules]

        # Create Pybel fingerprints for all the molecule bits. Store it in a tuple with the PK of the molecule
        all_fingerprints = [(pybel.Fingerprint(f[0]), f[1]) for f in all_bitlists]

        # Convert query molecue into a fingerprint
        query_molecule = pybel.readstring('smiles', query_smiles)
        query_fp = query_molecule.calcfp()

        # Compare each molecule fingerprint with the query fingerprint, store tuple if the score is above a threshold
        matches = [(query_fp | molecule[0], molecule[1]) for molecule in all_fingerprints
                   if query_fp | molecule[0] > tanimoto_threshold]

        ids = []
        if matches:
            scores, ids = zip(*matches)

        return self.filter(pk__in=ids)


class Molecule(Data):
    smiles = models.CharField(max_length=1000, verbose_name='SMILES', unique=True,
                              help_text="Benzene: C1=CC=CC=C1", validators=[validate_smiles])
    inchi = models.CharField(max_length=1000, verbose_name='InChI', unique=True)
    image = models.FileField(upload_to='molecules', verbose_name='Picture', blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)
    representation_3d = models.TextField(blank=True)

    fingerprint = models.CharField(max_length=1000, blank=True, null=True)

    objects = MoleculeManager()

    class Meta:
        unique_together = ('smiles', 'inchi')

    def __str__(self):
        return self.inchi


class Article(Data):
    author = models.CharField(max_length=1000)
    title = models.CharField(max_length=1000)
    journal = models.CharField(max_length=250)
    volume = models.CharField(max_length=100, blank=True, null=True)
    doi = models.CharField(max_length=500, verbose_name='DOI', unique=True)
    pages = models.CharField(max_length=20, blank=True, null=True)
    issue_nr = models.CharField(max_length=100, blank=True, null=True)
    eid = models.CharField(blank=True, null=True, max_length=100)
    year = models.DateField()
    electronic_id = models.CharField(max_length=250)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    keywords = models.CharField(max_length=1000, blank=True, null=True)

    def __str__(self):
        return '{} - "{}", vol. {}, issue {}, {} '.format(self.doi, self.title, self.volume, self.issue_nr, self.year)


class Spectrum(Data):
    absorption_maxima = models.DecimalField(blank=True, null=True, decimal_places=4, max_digits=10, help_text="[nm]")
    emission_maxima = models.DecimalField(blank=True, null=True, decimal_places=4, max_digits=10)
    solvent = models.CharField(max_length=100)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    # Foreignkey, not one to one
    molecule = models.OneToOneField(Molecule, related_name='spectrum')
    article = models.ForeignKey(Article, related_name='spectra')

    class Meta:
        unique_together = ('molecule', 'article')
        verbose_name = "Molecule spectrum"
        verbose_name_plural = "Molecule spectra"

    def __str__(self):
        return '{} - abs. max {}, emi. max {}'.format(self.molecule, self.absorption_maxima, self.emission_maxima)


class Performance(Data):
    voc = models.DecimalField(verbose_name='VOC', decimal_places=4, max_digits=15, help_text='[mV]')
    jsc = models.DecimalField(verbose_name='JSC', decimal_places=5, max_digits=15, help_text='[mA/cm^2]')
    ff = models.DecimalField(verbose_name='FF', decimal_places=5, max_digits=13)
    pce = models.DecimalField(verbose_name='PCE', decimal_places=5, max_digits=13, help_text='%, 0-1')
    electrolyte = models.CharField(max_length=1000)
    active_area = models.CharField(max_length=30, help_text='[cm2]', blank=True, null=True)
    co_adsorbent = models.CharField(max_length=250, blank=True, null=True)
    co_sensitizer = models.CharField(max_length=1000, blank=True, null=True)
    semiconductor = models.CharField(max_length=1000)
    dye_loading = models.CharField(max_length=1000, help_text='[nmol/cm2]', blank=True, null=True)
    exposure_time = models.CharField(max_length=500, blank=True, null=True)
    solar_simulator = models.CharField(max_length=1000, default='AM 1.5g')

    short_id = TinyUUIDField(length=10)

    # Now optional
    comment = models.TextField(blank=True, null=True)
    created = models.DateTimeField(auto_now_add=True, auto_now=False)

    article = models.ForeignKey(Article, related_name='performances')
    molecule = models.ForeignKey(Molecule, related_name='performances')
    contribution = GenericRelation(AtomicContribution)

    def __str__(self):
        return str(self.molecule)

    def get_absolute_url(self):
        return reverse("dye:performance-detail", kwargs={'short_id': self.short_id})

    class Meta:
        verbose_name = "DSSC performance"
        unique_together = ('article', 'molecule', 'voc', 'jsc', 'ff', 'pce')


class ContributionManager(models.Manager):
    def create_from_data(self, upload, *args, **kwargs):
        """
        Creates atomic uploads that can point to any Model: Article, spectra, molecule etc.
        A contribution is then assigned to all these atomic uploads
        """

        new_contribution = self.create(*args, **kwargs)
        unraveled_data_entries = list(chain.from_iterable(upload))
        content = [(e, ContentType.objects.get_for_model(e), created) for e, created in unraveled_data_entries]

        for entry_data, entry_type, created in content:
            if created:
                atom = AtomicContribution.objects.create(
                    content_type=entry_type,
                    object_id=entry_data.pk,
                )
                new_contribution.items.add(atom)
        return new_contribution


class Contribution(Data):
    user = models.ForeignKey(User, related_name='contributions')

    short_id = TinyUUIDField(length=10)
    created = models.DateTimeField(auto_now_add=True)
    status = models.PositiveSmallIntegerField(choices=APPROVAL_STATES, default=APPROVAL_STATES.WAITING)

    items = models.ManyToManyField(AtomicContribution, related_name='upload')

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
        for item in self.items.all():
            try:
                item.content_object.set_status(status)
                item.content_object.save()
            except AttributeError:
                pass
        self.status = status

    def save(self, *args, **kwargs):
        super(Contribution, self).save(*args, **kwargs)
        if self.pk:
            self.set_approval_state(self.status)

    def get_absolute_url(self):
        return reverse("dye:single-evaluation", kwargs={'short_id': self.short_id})

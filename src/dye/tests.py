from django.contrib.auth.models import User
from django.shortcuts import reverse
from django.test import TestCase, Client

from .models import Performance
from .validators import validate_inchi, validate_smiles, ValidationError

class FileUploadTestCase(TestCase):
    def setUp(self):
        self.username, self.email, self.password = 'testuser', 'test@test.test', 'testpassword'
        self.user = User.objects.create(username=self.username,
                                        email=self.email)
        self.user.set_password(self.password)
        self.user.save()
        self.client = Client()
        self.client.login(username=self.username, password=self.password)
        self.batch_url = reverse('dye:file-upload')
        self.single_upload_url = reverse('dye:single-upload')

    def test_missing_start_tag(self):
        with open('dye/testfiles/missing_start_tag.xlsx', 'rb') as fp:
            response = self.client.post(self.batch_url, {'file': fp}, follow=True)
            self.assertContains(response, 'Could not find start-tag. Compare your sheet with the sample sheet.')

    def test_missing_or_invalid_doi(self):
        with open('dye/testfiles/invalid_DOI.xlsx', 'rb') as fp:
            response = self.client.post(self.batch_url, {'file': fp}, follow=True)
            self.assertContains(response, 'DOI not found')

    def test_corrupt_file(self):
        with open('dye/testfiles/corrupt_file.notxlsx', 'rb') as fp:
            response = self.client.post(self.batch_url, {'file': fp}, follow=True)
            self.assertContains(response, 'The file was not recognized as a valid spreadsheet file. '
                                          'Please download the sample file and try again.')

    def test_missing_row_data(self):
        with open('dye/testfiles/missing_all_fields.xlsx', 'rb') as fp:
            response = self.client.post(self.batch_url, {'file': fp}, follow=True)
            self.assertContains(response, 'Critical error at row 8.')

    def test_invalid_performance_decimals(self):
        # Test that JSC, FF, VOC and PCE are valid numbers
        with open('dye/testfiles/invalid_decimal_fields.xlsx', 'rb') as fp:
            response = self.client.post(self.batch_url, {'file': fp}, follow=True)
            self.assertContains(response, 'Enter a number.', count=4)

    def test_pass_upload(self):
        with open('dye/testfiles/single_contribution.xlsx', 'rb') as fp:
            response = self.client.post(self.batch_url, {'file': fp}, follow=True)
            self.assertContains(response, 'The data was uploaded and is awaiting review. Thank you!')

    def test_single_upload(self):
        form_data = {
            # DOI
            'doi': '10.1039/c5ra02720a',

            # Performance
            'voc': '700', 'jsc': '15.00', 'ff': '0.5', 'pce': '3.14', 'electrolyte': 'test_electrolyte',
            'active_area': 'test_active_area', 'co_adsorbent': 'test_co_adsorbent',
            'co_sensitizer': 'test_co_sensitizer',
            'semiconductor': 'test_semiconductor', 'dye_loading': 'test_dye_loading',
            'exposure_time': 'test_exposure_time',
            'solar_simulator': 'test_solar_simulator', 'keywords': 'test_keywords', 'comment': 'test_comment',

            # Spectrum
            'absorption_maxima': '1', 'emission_maxima': '1',
            'solvent': 'test_solvent',

            # Molecule
            'smiles': 'O=Cc1ccc(O)c(OC)c1', 'inchi': 'MWOOGOJBHIARFG-UHFFFAOYSA-N'
        }
        response = self.client.post(self.single_upload_url, form_data, follow=True)
        self.assertContains(response, 'The data was uploaded and is awaiting review. Thank you!')


class SingleUploadTest(TestCase):
    def setUp(self):
        self.username, self.email, self.password = 'testuser', 'test@test.test', 'testpassword'
        self.user = User.objects.create(username=self.username,
                                        email=self.email)
        self.user.set_password(self.password)
        self.user.save()
        self.client = Client()
        self.client.login(username=self.username, password=self.password)
        self.single_upload_url = reverse('dye:single-upload')

    def test_data_integrity(self):
        form_data = {
            # DOI
            'doi': '10.1039/c5ra02720a',

            # Performance
            'voc': '700', 'jsc': '15.00', 'ff': '0.5', 'pce': '3.14', 'electrolyte': 'test_electrolyte',
            'active_area': 'test_active_area', 'co_adsorbent': 'test_co_adsorbent',
            'co_sensitizer': 'test_co_sensitizer',
            'semiconductor': 'test_semiconductor', 'dye_loading': 'test_dye_loading',
            'exposure_time': 'test_exposure_time',
            'solar_simulator': 'test_solar_simulator', 'keywords': 'test_keywords', 'comment': 'test_comment',

            # Spectrum
            'absorption_maxima': '1', 'emission_maxima': '1',
            'solvent': 'test_solvent',

            # Molecule
            'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'inchi': 'RYYVLZVUVIJVGH-UHFFFAOYSA-N'
        }
        response = self.client.post(self.single_upload_url, form_data, follow=True)
        self.assertContains(response, 'The data was uploaded and is awaiting review. Thank you!')

        performance = Performance.objects.last()
        self.assertEquals(performance.molecule.smiles, form_data.get('smiles'))

        self.assertIn('OpenBabel', performance.molecule.representation_3d)


class ValidationTest(TestCase):
    def test_invalid_inchi(self):
        with self.assertRaises(ValidationError):
            validate_smiles('E')

    def test_valid_inchi(self):
        self.assertIsNone(validate_smiles('C'))

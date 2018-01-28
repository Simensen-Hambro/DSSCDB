from django.core.urlresolvers import reverse
from django.test import Client
from django.test import TestCase

from usermanagement.models import *
import os

class UserTokenTestCase(TestCase):
    def setUp(self):
        # Creates a user that is not yet activated
        self.user = User.objects.create(username='test')
        self.user.is_active = False
        self.user.save()

        # Cretes three different UserToken objects used in other tests
        self.expired_token = UserToken.objects.create(user=self.user)
        expired_date = timezone.now() - timedelta(hours=VALID_TIME + 1)
        self.expired_token.created = expired_date
        self.expired_token.save()

        self.activation_token = UserToken.objects.create(user=self.user)
        self.set_password_token = UserToken.objects.create(user=self.user)

    def test_expiration_of_expired_token(self):
        self.assertTrue(self.expired_token.expired(), 'Checks if UserToken expires')

    def test_activation_of_user(self):
        self.assertFalse(self.user.is_active, 'Checks if user is active before activation')
        self.activation_token.activate()
        self.assertTrue(self.user.is_active, 'Checks if user is active before activation')

        # Checks if the UserToken does exist after activation
        with self.assertRaises(UserToken.DoesNotExist):
            UserToken.objects.get(key=self.activation_token.key)

    def test_setting_password_with_token(self):
        password = 'test_password'
        self.assertFalse(self.user.check_password(password), 'Checks if password matches before its been set')
        self.set_password_token.set_password(password)
        self.assertTrue(self.user.check_password(password), 'Checks if password matches after its been set')

        # Checks if the UserToken does exist after activation
        with self.assertRaises(UserToken.DoesNotExist):
            UserToken.objects.get(key=self.set_password_token.key)

    def test_prune_expired_tokens(self):
        # Check is UserToken exists before deletion
        self.assertEqual(self.expired_token, UserToken.objects.get(key=self.expired_token.key))

        # Deletes all the expired tokens
        UserToken.objects.prune_expired()

        # Checks that the expired UserToken does not exist after deletion
        with self.assertRaises(UserToken.DoesNotExist):
            UserToken.objects.get(key=self.expired_token.key)


class LoginViewAndLogoutViewTestCase(TestCase):
    def setUp(self):
        self.user = User.objects.create(username='test', email='test@test.test')
        self.password = 'Test1234'
        self.user.set_password(self.password)
        self.user.save()
        self.url = reverse('user:login')

        # Testing with ReCaptcha requires a setting and passing a flag to the form
        # https://github.com/praekelt/django-recaptcha#unit-testing
        os.environ['RECAPTCHA_TESTING'] = 'True'
        self.payload = {
            'email': self.user.email,
            'password': self.password,
            'g-recaptcha-response': 'PASSED'
        }
        self.client = Client()

    def tearDown(self):
        os.environ['RECAPTCHA_TESTING'] = 'False'

    def test_login_succesful(self):
        self.client.get(self.url)
        response = self.client.post(self.url, self.payload)
        self.assertRedirects(response, reverse('index'))

    def test_login_while_already_logged_in(self):
        self.client.post(self.url, self.payload)
        response = self.client.post(self.url, self.payload)
        self.assertEqual(response.status_code, 404, 'Should get 404 when trying to login while already logged in')

    def test_login_with_next_paramenter(self):
        self.client.logout()
        payload = self.payload
        payload['next'] = reverse('index')
        response = self.client.post(self.url, payload)
        self.assertRedirects(response, reverse('index'))

    def test_login_with_not_existing_user(self):
        self.client.logout()
        payload = self.payload
        payload['email'] = 'notExistingEmail@test.test'
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])

    def test_login_with_inactive_user(self):
        self.user.is_active = False
        self.user.save()
        response = self.client.post(self.url, self.payload)
        self.assertIsNotNone(response.context['errors'])
        self.user.is_active = True
        self.user.save()

    def test_login_with_wrong_password(self):
        payload = self.payload
        payload['password'] = 'WrongPassword123'
        response = self.client.post(self.url, self.payload)
        self.assertIsNotNone(response.context['errors'])

    def test_logout(self):
        response = self.client.post(reverse('user:logout'))
        self.assertRedirects(response, reverse('index'))


class SignupUserViewTestCase(TestCase):
    fixtures = ['test_fixtures.json', ]

    def setUp(self):
        self.username = 'test'
        self.email = 'test@test.test'
        self.password = 'Test1234'
        self.url = reverse('user:signup')
        os.environ['RECAPTCHA_TESTING'] = 'True'
        self.payload = {
            'first_name': 'test',
            'last_name': 'test',
            'email': self.email,
            'organization': 'NTNU',
            'affiliation': 'Professor',
            'password': self.password,
            'confirm_password': self.password,
            'g-recaptcha-response': 'PASSED',
            'agree_to_terms': "True",
        }
        self.client = Client()

    def tearDown(self):
        os.environ['RECAPTCHA_TESTING'] = 'False'

    def test_successful_signup(self):
        self.client.get(self.url)
        response = self.client.post(self.url, self.payload)
        self.assertRedirects(response, reverse('index'))

    def test_signup_with_existing_email(self):
        self.client.post(self.url, self.payload)
        response = self.client.post(self.url, self.payload)
        self.assertIsNotNone(response.context['errors'])

    def test_signup_with_non_matching_passwords(self):
        payload = self.payload
        payload['confirm_password'] = self.password + "someMore"
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])

    def test_signup_with_invalid_password(self):
        payload = self.payload
        payload['password'] = "passwordWithoutNumber"
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])


class ForgotPasswordViewTestCase(TestCase):
    fixtures = ['test_fixtures.json', ]

    def setUp(self):
        self.user = User.objects.create(username='test', email="test@test.test")
        self.user.set_password("Test1234")
        self.user.save()
        self.url = reverse('user:forgot-password')

        os.environ['RECAPTCHA_TESTING'] = 'True'
        self.payload = {
            'email': self.user.email,
            'g-recaptcha-response': 'PASSED'
        }
        self.client = Client()


    def tearDown(self):
        os.environ['RECAPTCHA_TESTING'] = 'False'

    def test_forgot_password_successful(self):
        self.client.get(self.url)
        response = self.client.post(self.url, self.payload)
        self.assertRedirects(response, reverse('index'))

    def test_forgot_password_on_inactive_user(self):
        self.user.is_active = False
        self.user.save()
        response = self.client.post(self.url, self.payload)
        self.assertRedirects(response, reverse('index'))
        self.user.is_active = True
        self.user.save()

    def test_forgot_password_with_non_existing_email(self):
        payload = self.payload
        payload['email'] = "nonExistingEmail@test.test"
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])


class ChangePasswordViewTestCase(TestCase):
    def setUp(self):
        self.username = 'test'
        self.user = User.objects.create(username=self.username,
                                        email='test@test.test')
        self.password = 'test'
        self.new_password = 'Test1234'
        self.user.set_password(self.password)
        self.user.save()
        self.url = reverse('user:change-password')
        self.payload = {
            'current_password': self.password,
            'new_password': self.new_password,
            'confirm_password': self.new_password,
        }

        self.client = Client()

    def test_change_password_without_authentication(self):
        response = self.client.post(self.url)
        self.assertEqual(response.status_code, 404,
                         'Should get 404 when trying to change password without being logged in')

    def test_change_password_successful(self):
        self.client.login(username=self.username, password=self.password)
        self.client.get(self.url)
        response = self.client.post(self.url, self.payload)
        self.assertRedirects(response, reverse('index'))

    def test_change_password_with_wrong_existing_password(self):
        self.client.login(username=self.username, password=self.password)
        payload = self.payload
        payload['current_password'] = 'wrong_password123'
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])

    def test_change_password_with_non_matching_passwords(self):
        self.client.login(username=self.username, password=self.password)
        payload = self.payload
        payload['confirm_password'] = 'wrong_password123'
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])

    def test_change_password_with_invalid_password(self):
        self.client.login(username=self.username, password=self.password)
        payload = self.payload
        payload['new_password'] = 'wrong'
        payload['confirm_password'] = 'wrong'
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])


class SetPasswordViewTestCase(TestCase):
    def setUp(self):
        self.username = 'test'
        self.user = User.objects.create(username=self.username,
                                        email='test@test.test')
        self.password = 'test'
        self.new_password = 'Test1234'
        self.url = reverse('user:set-password', args=[UserToken.objects.create(user=self.user).key])
        self.user.set_password(self.password)
        self.user.save()
        self.payload = {
            'new_password': self.new_password,
            'confirm_password': self.new_password,
        }

        self.client = Client()

    def test_set_password_without_token(self):
        invalid_url = reverse('user:set-password', args=["12345678123456781234567812345678"])  # Token as argument
        response = self.client.post(invalid_url)
        self.assertEqual(response.status_code, 404,
                         'Should get 404 when trying to set password without valid token')

    def test_set_password_successful(self):
        self.client.get(self.url)
        response = self.client.post(self.url, self.payload)
        self.assertRedirects(response, reverse('index'))

    def test_set_invalid_password(self):
        payload = self.payload
        payload['new_password'] = 'wrong'
        response = self.client.post(self.url, payload)
        self.assertIsNotNone(response.context['errors'])


class ActivateUserTestCase(TestCase):
    fixtures = ['email_templates.json']
    def setUp(self):
        self.username = 'test'
        self.user = User.objects.create(username=self.username,
                                        email='test@test.test')
        self.url = reverse('user:activate', args=[UserToken.objects.create(user=self.user).key])
        self.client = Client()

    def test_activation_without_valid_token(self):
        invalid_url = reverse('user:activate', args=["12345678123456781234567812345678"])
        response = self.client.post(invalid_url)
        self.assertEqual(response.status_code, 404,
                         'Should get 404 when trying to activate without valid token')

    def test_valid_activation(self):
        response = self.client.post(self.url)
        self.assertRedirects(response, reverse('index'))

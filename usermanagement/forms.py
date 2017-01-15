from django.contrib.auth.models import User
from django.contrib.auth import login
from django import forms
from usermanagement.models import Profile


class LoginForm(forms.Form):
    email = forms.EmailField(label='Email', max_length=100)
    password = forms.CharField(label='Password', max_length=100, widget=forms.PasswordInput)

    def clean(self):
        cleaned_data = super(LoginForm, self).clean()
        email = cleaned_data.get("email")
        password = cleaned_data.get("password")

        error_msg = "Email or password is wrong"
        if email and password:
            try:
                user = User.objects.get(email=email)
                if not user.is_active:
                    self.add_error('email',
                                   'The account is not confirmed. Click forgot password to receive a new confirmation email')
                    return
            except User.DoesNotExist:
                self.add_error('email', error_msg)
                self.add_error('password', error_msg)
                return

            if not user.check_password(password):
                self.add_error('email', error_msg)
                self.add_error('password', error_msg)

    def authenticate_user(self, request):
        email = self.cleaned_data.get("email")
        try:
            user = User.objects.get(email=email)
            login(request, user)
        except User.DoesNotExist:
            print("Could not login user with email: {}".format(email))


class SignUpForm(forms.Form):
    first_name = forms.CharField(label='First name', max_length=255)
    last_name = forms.CharField(label='Last name', max_length=255)
    email = forms.EmailField(label='Email', max_length=255)
    organization = forms.CharField(label='Organization', max_length=255)
    affiliation = forms.CharField(label='Affiliation', max_length=255)
    password = forms.CharField(label='Password', max_length=255, widget=forms.PasswordInput)
    confirm_password = forms.CharField(label='Confirm password', max_length=255, widget=forms.PasswordInput)

    def clean(self):
        cleaned_data = super(SignUpForm, self).clean()
        email = cleaned_data.get("email")
        password = cleaned_data.get('password')
        confirm_password = cleaned_data.get('confirm_password')

        try:
            User.objects.get(email=email)
            self.add_error('email', 'The email address is already registered')
        except User.DoesNotExist:
            pass

        if password != confirm_password:
            msg = 'Passwords does not match'
            self.add_error('password', msg)
            self.add_error('confirm_password', msg)

        if not check_password(password):
            self.add_error('password',
                           'Password must be at least 8 characters, contain at least one number and at least one capital letter')

    def create_user(self):
        first_name = self.cleaned_data.get('first_name')
        last_name = self.cleaned_data.get('last_name')
        email = self.cleaned_data.get('email')
        organization = self.cleaned_data.get('organization')
        affiliation = self.cleaned_data.get('affiliation')
        password = self.cleaned_data.get('password')
        user = User.objects.create(username=email,
                                   email=email,
                                   first_name=first_name,
                                   last_name=last_name,
                                   is_active=False)
        user.set_password(password)
        user.save()

        Profile.objects.create(user=user,
                               organization=organization,
                               affiliation=affiliation)


class ForgotPasswordForm(forms.Form):
    email = forms.CharField(label='Email', max_length=100)

    def clean(self):
        cleaned_data = super(ForgotPasswordForm, self).clean()
        email = cleaned_data.get('email')

        try:
            user = User.objects.get(email=email)
        except User.DoesNotExist:
            self.add_error('email', 'The email is not registered')


class SetPasswordForm(forms.Form):
    new_password = forms.CharField(label='New password', max_length=100, widget=forms.PasswordInput)
    confirm_password = forms.CharField(label='Confirm new password', max_length=100, widget=forms.PasswordInput)

    def clean(self):
        cleaned_data = super(SetPasswordForm, self).clean()
        new_password = cleaned_data.get('new_password')
        confirm_password = cleaned_data.get('confirm_password')

        if new_password != confirm_password:
            self.add_error('new_password', 'Passwords does not match')
            self.add_error('confirm_password', 'Passwords does not match')

        if not check_password(new_password):
            self.add_error('new_password',
                           'Password must be at least 8 characters, contain at least one number and at least one capital letter')


class ChangePasswordForm(forms.Form):
    current_password = forms.CharField(label='Current password', max_length=100, widget=forms.PasswordInput)
    new_password = forms.CharField(label='New password', max_length=100, widget=forms.PasswordInput)
    confirm_password = forms.CharField(label='Confirm new password', max_length=100, widget=forms.PasswordInput)

    def __init__(self, email, *args, **kwargs):
        self.email = email
        super(ChangePasswordForm, self).__init__(*args, **kwargs)

    def clean(self):
        cleaned_data = super(ChangePasswordForm, self).clean()
        current_password = cleaned_data.get('current_password')
        new_password = cleaned_data.get('new_password')
        confirm_password = cleaned_data.get('confirm_password')
        email = self.email

        user = User.objects.get(email=email)

        if not user.check_password(current_password):
            self.add_error('current_password', 'Wrong password')

        if new_password != confirm_password:
            self.add_error('new_password', 'Passwords does not match')
            self.add_error('confirm_password', 'Passwords does not match')

        if not check_password(new_password):
            self.add_error('new_password',
                           'Password must be at least 8 characters, contain at least one number and at least one capital letter')

    def change_password(self, request):
        email = self.email
        new_password = self.cleaned_data.get('new_password')
        user = User.objects.get(email=email)
        user.set_password(new_password)
        user.save()
        login(request, user)


def check_password(password):
    # To short password
    if len(password) < 8:
        return False

    # No capital letter
    if password.lower() == password:
        return False

    # Check for at least a digit
    digit = False
    for ch in password:
        if ch.isdigit():
            digit = True

    if not digit:
        return False

    return True

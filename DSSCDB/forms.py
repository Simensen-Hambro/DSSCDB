from django import forms
from captcha.fields import ReCaptchaField

class ContactForm(forms.Form):
    contact_name = forms.CharField(required=False, label="Name")
    contact_email = forms.EmailField(required=False, label="Email")
    content = forms.CharField(required=True, widget=forms.Textarea, label="Message")
    captcha = ReCaptchaField()

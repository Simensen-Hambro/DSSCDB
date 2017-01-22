from .forms import ContactForm
from django.contrib import messages
from django.conf import settings
from django.shortcuts import reverse, redirect, render
from post_office import mail


def index(request):
    context = {}

    return render(request, 'DSSCDB/index.html', context)


def contact_us(request):
    form_data = ContactForm(request.POST or None)

    if form_data.is_valid():
        messages.add_message(request, messages.SUCCESS, 'The message has been received. Thanks for contacting us!',
                             extra_tags="Received!")

        _, mail_to = zip(*settings.ADMINS)

        mail.send(
            ['carljh@stud.ntnu.no'],
            subject='Contact form submission',
            message='You have been contacted! '
                    'Who: {{ contact_name|default_if_none:"Not given" }} '
                    'Email: {{ contact_email|default_if_none:"Not given" }} '
                    'When: {% now "j. F - H:i" %} '
                    'Message: {{message}}',
            context={'message': form_data.cleaned_data.get('content'),
                     'contact_name': form_data.cleaned_data.get('contact_name'),
                     'contact_email': form_data.cleaned_data.get('contact_email')
                     },
        )
        return redirect(reverse('index'))
    else:
        return render(request, 'DSSCDB/contact-us.html', {
            'contact_form': form_data,
        })

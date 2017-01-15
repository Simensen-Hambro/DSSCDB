from django.http import HttpResponseRedirect
from django.contrib.auth import logout
from usermanagement.forms import *
from django.shortcuts import reverse, render, Http404
from django.contrib import messages
from usermanagement.models import *
from post_office import mail
from django.conf import settings


def login_view(request):
    if request.method == 'POST':
        form = LoginForm(request.POST)
        if form.is_valid():
            form.authenticate_user(request)
            return HttpResponseRedirect(reverse('index'))
    else:
        form = LoginForm()

    context = {
        'form': form
    }

    return render(request, 'usermanagement/login.html', context)


def logout_view(request):
    logout(request)
    return HttpResponseRedirect(reverse('index'))


def signup_view(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            form.create_user()
            user = User.objects.get(email=form.cleaned_data.get('email'))
            token = UserToken.objects.create(user=user)
            mail.send(recipients=[user.email],
                      sender=settings.DEFAULT_FROM_MAIL,
                      template='activation_email',
                      context={'request': request, 'user': user, 'token': token})
            messages.add_message(request, messages.INFO, 'You will receive a confirmation email to verify you email address')
            return HttpResponseRedirect(reverse('index'))

    else:
        form = SignUpForm()

    context = {
        'form': form
    }

    return render(request, 'usermanagement/signup.html', context)


def forgot_password_view(request):
    if request.method == 'POST':
        form = ForgotPasswordForm(request.POST)
        if form.is_valid():
            user = User.objects.get(email=form.cleaned_data.get('email'))
            token = UserToken.objects.create(user=user)
            if not user.is_active:
                mail.send(recipients=[user.email],
                          sender=settings.DEFAULT_FROM_MAIL,
                          template='activation_email',
                          context={'request': request, 'user': user, 'token': token})
                messages.add_message(request, messages.INFO,
                                     'You will receive a confirmation email to verify you email address')
            else:
                mail.send(recipients=[user.email],
                          sender=settings.DEFAULT_FROM_MAIL,
                          template='set_password_email',
                          context={'request': request, 'user': user, 'token': token})
                messages.add_message(request, messages.INFO,
                                     'You will receive an email with further instruction to reset your password')

            return HttpResponseRedirect(reverse('index'))
    else:
        form = ForgotPasswordForm()

    context = {
        'form': form
    }

    return render(request, 'usermanagement/forgot_password.html', context)


def change_password_view(request):
    if not request.user.is_authenticated:
        raise Http404

    if request.method == 'POST':
        form = ChangePasswordForm(data=request.POST, email=request.user.email)
        if form.is_valid():
            form.change_password(request)
            messages.add_message(request, messages.INFO, 'The password is now changed')
            return HttpResponseRedirect(reverse('index'))

    else:
        form = ChangePasswordForm(email=request.user.email)

    context = {
        'form': form
    }
    return render(request, 'usermanagement/change_password.html', context)


def activate(request, key):
    try:
        token = UserToken.objects.get(key=key)
        token.activate()
        messages.add_message(request, messages.INFO, 'Your account is now activated')
        return HttpResponseRedirect(reverse('index'))
    except UserToken.DoesNotExist:
        raise Http404


def set_password_view(request, key):
    try:
        token = UserToken.objects.get(key=key)
    except UserToken.DoesNotExist:
        raise Http404

    if request.method == 'POST':
        form = SetPasswordForm(request.POST)
        if form.is_valid():
            token.set_password(form.cleaned_data.get('new_password'))
            messages.add_message(request, messages.INFO, 'The new password is now set')
            return HttpResponseRedirect(reverse('index'))

    else:
        form = SetPasswordForm()

    context = {
        'form': form,
        'key': key,
    }
    return render(request, 'usermanagement/set_password.html', context)


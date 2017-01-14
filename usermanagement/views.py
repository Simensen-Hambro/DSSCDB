from django.http import HttpResponseRedirect
from django.contrib.auth import logout
from usermanagement.forms import *
from django.shortcuts import reverse, render, Http404
from django.contrib import messages


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

            messages.add_message(request, messages.SUCCESS, 'The password is now changed')
            return HttpResponseRedirect(reverse('index'))

    else:
        form = ChangePasswordForm(email=request.user.email)

    context = {
        'form': form
    }
    return render(request, 'usermanagement/change_password.html', context)


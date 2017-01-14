from django.http import HttpResponseRedirect
from django.contrib.auth import logout
from usermanagement.forms import LoginForm, SignUpForm
from django.shortcuts import reverse, render


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

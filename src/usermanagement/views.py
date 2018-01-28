from django.conf import settings
from django.contrib import messages
from django.contrib.auth import logout
from django.contrib.auth.decorators import permission_required
from django.http import HttpResponseRedirect
from django.shortcuts import reverse, render, Http404
from post_office import mail

from usermanagement.forms import *
from usermanagement.models import *
from .helpers import get_ip_address_from_request


def login_view(request):
    if request.user.is_authenticated:
        raise Http404
    if request.method == 'POST':
        next = request.POST.get('next')
        form = LoginForm(request.POST)
        if form.is_valid():
            form.authenticate_user(request)
            if request.POST.get('next') == None or request.POST.get('next') == 'None':
                return HttpResponseRedirect(reverse('index'))
            else:
                return HttpResponseRedirect(request.POST.get('next'))
    else:
        form = LoginForm()
        next = request.GET.get('next')

    context = {
        'form': form,
        'next': next
    }

    return render(request, 'usermanagement/login.html', context)


def logout_view(request):
    logout(request)
    return HttpResponseRedirect(reverse('index'))


def signup_view(request):
    if request.method == 'POST':
        form = SignUpForm(request.POST)
        if form.is_valid():
            ip = get_ip_address_from_request(request)
            form.create_user(ip)
            user = User.objects.get(email=form.cleaned_data.get('email'))
            token = UserToken.objects.create(user=user)
            mail.send(recipients=[user.email],
                      sender=settings.DEFAULT_FROM_MAIL,
                      template='activation_email',
                      context={'request': request, 'user': user, 'token': token})

            messages.add_message(request, messages.INFO,
                                 'You will receive a confirmation email to verify you email address. '
                                 'An administrator will then approve your account before it can be used.')
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
        messages.add_message(request, messages.INFO, 'Your email address is now confirmed. '
                                                     'An administrator will approve your account before it can be used.')
        approval_token = UserApprovalToken.objects.create(user=token.user)
        _, mail_to = zip(*settings.ADMINS)
        mail.send(recipients=list(mail_to),
                  sender=settings.DEFAULT_FROM_MAIL,
                  template='user_registration_notification',
                  context={'user': token.user, 'token': approval_token, 'request': request})
        return HttpResponseRedirect(reverse('index'))
    except UserToken.DoesNotExist:
        raise Http404


@permission_required('usermanagement.change_user')
def approve_user(request, key):
    try:
        token = UserApprovalToken.objects.get(key=key)
        token.activate()
        user = token.user
        messages.add_message(request, messages.INFO, 'User has been approved')
        mail.send(recipients=[user.email],
                  sender=settings.DEFAULT_FROM_MAIL,
                  template='user_account_approval',
                  context={'request': request, 'user': user, 'token': token})
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


def profile(request):
    if not request.user.is_authenticated:
        raise Http404

    context = {
        'user_profile': request.user.profile,
    }

    return render(request, 'usermanagement/profile.html', context)


def admin_users(request):
    if not request.user.is_staff:
        raise Http404

    context = {
        'profiles': Profile.objects.all().order_by('organization', 'user__first_name', 'user__last_name'),
    }

    return render(request, 'usermanagement/admin_users.html', context)

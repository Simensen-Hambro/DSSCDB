from django.db import models
from django.contrib.auth.models import User
from extended_choices import Choices
from django.utils import timezone
from datetime import timedelta
import uuid
from django.contrib.auth.models import Permission

# Time the activation is valid, in hours
VALID_TIME = 48


class Profile(models.Model):
    user = models.OneToOneField(User, related_name='profile')

    organization = models.CharField(max_length=100)
    affiliation = models.CharField(max_length=250)

    def __str__(self):
        return '{} {}'.format(self.user.first_name, self.user.last_name)


class RequestUploadUserStatus(models.Model):
    STATES = Choices(
        ('WAITING', 1, 'Waiting'),
        ('APPROVED', 2, 'Approved'),
        ('DENIED', 3, 'Denied'),
    )
    user = models.ForeignKey(User)
    status = models.PositiveSmallIntegerField(choices=STATES, default=STATES.WAITING)
    created = models.DateTimeField(auto_now_add=True)
    edited = models.DateTimeField(auto_now=True)

    def approve(self):
        upload_permission = Permission.objects.get(codename='upload_performance_data')
        self.status = self.STATES.APPROVED
        self.user.user_permissions.add(upload_permission)

    def deny(self):
        self.status = self.STATES.DENY
        self.user.user_permissions.clear()

    def __str__(self):
        return '{} {} - {}'.format(self.user.first_name, self.user.last_name,
                                   self.STATES.for_value(self.status).display)

    def save(self, *args, **kwargs):
        # Whenever the object is saved, set the new permissions
        if self.status == self.STATES.APPROVED:
            self.approve()
        elif self.status == self.STATES.DENIED:
            self.deny()
        super(RequestUploadUserStatus, self).save()

    class Meta:
        verbose_name = "Upload status request"
        verbose_name_plural = "Upload status requests"


class UserTokenManager(models.Manager):
    def prune_expired(self):
        self.filter(created__lt=timezone.now() - timedelta(hours=VALID_TIME)).delete()


class UserToken(models.Model):
    user = models.ForeignKey(User)
    key = models.UUIDField(default=uuid.uuid4, editable=False)
    created = models.DateTimeField(auto_now=False, auto_now_add=True)

    objects = UserTokenManager()

    # Activates the user and deletes the authentication object
    def activate(self):
        self.user.is_active = True
        self.user.save()
        self.delete()

    # Set the password and deletes the authentication object
    def set_password(self, password):
        self.user.set_password(password)
        self.user.save()
        self.delete()

    # Checks if the authentication object is expired
    def expired(self):
        return not timezone.now() < timedelta(hours=VALID_TIME) + self.created

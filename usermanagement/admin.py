from django.contrib import admin

from .models import Profile, RequestUploadUserStatus

admin.site.register(Profile)
admin.site.register(RequestUploadUserStatus)
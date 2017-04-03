from django.contrib import admin

from .models import Profile, RequestUploadUserStatus, UserToken

admin.site.register(Profile)
admin.site.register(RequestUploadUserStatus)
admin.site.register(UserToken)

from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^upload', views.upload_data, name='upload'),
]

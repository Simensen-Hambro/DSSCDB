from django.conf.urls import url

from dye import views

urlpatterns = [
    url(r'^upload_file', views.file_upload, name='file-upload'),
    url(r'^upload', views.single_upload, name='single-upload'),
    url(r'^list-performance', views.list_performance, name='list-performance'),
]

from django.conf.urls import url

from dye import views

urlpatterns = [
    url(r'^upload_file', views.file_upload, name='file-upload'),
    url(r'^upload', views.single_upload, name='single-upload'),
    url(r'^performance-list', views.performance_list, name='performance-list'),
    url(r'^performance-detail/(?P<short_id>[0-9a-f]{8})/$', views.performance_details, name='performance-detail'),
]

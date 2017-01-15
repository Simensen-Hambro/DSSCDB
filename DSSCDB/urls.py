from django.conf.urls import url, include
from django.contrib import admin
from DSSCDB import views
from usermanagement import views as usermanagement_views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^admin/', admin.site.urls),
    url(r'^user/', include('usermanagement.urls', namespace='user')),
    url(r'^hub/', include('dye.urls', namespace='dye')),
]

from django.conf.urls import url, include
from django.contrib import admin
from DSSCDB import views
from django.conf import settings
from django.conf.urls.static import static
from usermanagement import views as u_views
from django.contrib.staticfiles import views as static_views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^admin/login/', u_views.login_view),
    url(r'^admin/', admin.site.urls, ),
    url(r'^user/', include('usermanagement.urls', namespace='user')),
    url(r'^hub/', include('dye.urls', namespace='dye')),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

urlpatterns += [
    url(r'^', include('django.contrib.flatpages.urls', namespace='flatpages')),
]

from django.conf.urls import url
from usermanagement import views

urlpatterns = [
    url(r'^login/$', views.login_view, name='login'),
    url(r'^logout/$', views.logout_view, name='logout'),
    url(r'^signup/$', views.signup_view, name='signup'),
    url(r'^forgot-password/$', views.forgot_password_view, name='forgot-password'),
    url(r'^set-password/(?P<key>[0-9a-z-]+)/$', views.set_password_view, name='set-password'),
    url(r'^change-password/$', views.change_password_view, name='change-password'),
    url(r'^activate/(?P<key>[0-9a-z-]+)/$', views.activate, name='activate'),
    url(r'^approve/(?P<key>[0-9a-z-]+)/$', views.approve_user, name='approve'),
    url(r'^profile/$', views.profile, name='profile'),
    url(r'^admin-users/$', views.admin_users, name='admin-users'),
]

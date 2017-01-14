from django.conf.urls import url
from usermanagement import views

urlpatterns = [
    url(r'^login/$', views.login_view, name='login'),
    url(r'^logout/$', views.logout_view, name='logout'),
    url(r'^signup/$', views.signup_view, name='signup'),
    url(r'^forgot-password/$', views.forgot_password_view, name='forgot-password'),
    url(r'^change-password/$', views.change_password_view, name='change-password'),
]

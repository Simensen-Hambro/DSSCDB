import os

from django.core.management import call_command
from django.core.wsgi import get_wsgi_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "DSSCDB.settings")

application = get_wsgi_application()
from django.conf import settings

try:
    import uwsgidecorators

    @uwsgidecorators.timer(60)
    def send_queued_mail(num):
        """Send queued mail every 60 seconds"""
        call_command('send_queued_mail', processes=1)


    @uwsgidecorators.cron(4, 4, -1, -1, 1)
    def send_database_backup_email(num):
        call_command('send_database_backup', **{'--backup_password': settings.DATABASE_EMAIL_PASSWORD})

except ImportError:
    print("uwsgidecorators not found. Cron and timers are disabled")

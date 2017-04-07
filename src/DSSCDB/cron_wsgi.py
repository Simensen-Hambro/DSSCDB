import os

from django.core.management import call_command
from django.core.wsgi import get_wsgi_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "DSSCDB.settings")

application = get_wsgi_application()

try:
    import uwsgidecorators
    

    @uwsgidecorators.timer(10)
    def send_queued_mail(num):
        """Send queued mail every 10 seconds"""
        call_command('send_queued_mail', processes=1)

except ImportError:
    print("uwsgidecorators not found. Cron and timers are disabled")

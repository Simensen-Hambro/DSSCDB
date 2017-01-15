from post_office import mail
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    args = ''
    help = 'Test send mail to an address'

    def handle(self, *args, **options):
        send_mail(args[0])


def send_mail(email):
    mail.send(recipients=[email],
              sender='carl.j.v.hambro@ntnu.no',
              subject='Test',
              message='Melding')

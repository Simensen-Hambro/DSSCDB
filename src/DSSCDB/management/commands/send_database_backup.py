import os
import time
from subprocess import Popen, PIPE

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from post_office.mail import send
from post_office.models import EmailTemplate


class Command(BaseCommand):
    help = "Dumps the database to a encrypted archive before it is sent to the administrators"

    def add_arguments(self, parser):
        parser.add_argument('--destination', dest='backup_directory', default='backups',
                            help='Destination (path) where to place database dump file.')
        parser.add_argument('--quiet', dest='quiet', action='store_true', default=False, help='Be silent.')
        parser.add_argument('--debug', dest='debug', action='store_true', default=False, help='Be silent.')
        parser.add_argument('--db-name', dest='database_name', default='default',
                            help='Name of database (as defined in settings.DATABASES[]) to dump.')
        parser.add_argument('--backup_password', dest='backup_password', default='change_me',
                            help='Password for the backup file')

    def handle(self, *args, **options):
        db_name = options.get('database_name') or settings.DATABASES.get('default')
        if not db_name:
            raise CommandError('Default database not found, please specify. ')

        database = settings.DATABASES.get(db_name)

        debug = options.get('debug')
        self.engine = database.get('ENGINE')
        self.db = database.get('NAME', '') or ''
        self.user = database.get('USER') or ''
        self.password = database.get('PASSWORD') or ''
        self.host = database.get('HOST')
        self.port = database.get('PORT')
        self.quiet = database.get('quiet')
        self.backup_password = options.get('backup_password')

        if not ('postgres' in self.engine):
            raise CommandError('Only Postgres database is supported.')

        backup_directory = options.get('backup_directory')

        outfile = self.destination_filename(backup_directory, self.db)
        command = self.get_postgresql_command(outfile)
        dumpfile = Popen("PGPASSWORD={} {}".format(self.password, command), shell=True, stdout=PIPE, stderr=PIPE,
                         executable='/bin/bash')
        dumpfile.wait()
        zip_file = outfile + '.zip'
        zip_command = 'zip --password {} {} {}'.format(self.backup_password, zip_file, outfile)
        encryptfile = Popen(zip_command,
                            shell=True, stdout=PIPE, stderr=PIPE, executable='/bin/bash')
        encryptfile.wait()
        self.send_email(zip_file)
        if debug:
            self.stdout.write(self.style.NOTICE('{}\n{}'.format(command, zip_command)))

        self.stdout.write(self.style.SUCCESS('Backup sent!'))

    def destination_filename(self, backup_directory, database_name):
        return os.path.join(backup_directory,
                            'backup_{}_backup_{}.sql'.format(database_name, time.strftime('%Y%m%d-%H%M%S')))

    def get_postgresql_command(self, outfile):
        if not self.quiet:
            print('Doing PostgreSQL backup of database "{}" into {}'.format(self.db, outfile))

        main_args = []
        if self.user:
            main_args += ['--username={}'.format(self.user)]
        if self.host:
            main_args += ['--host={}'.format(self.host)]
        if self.port:
            main_args += ['--port={}'.format(self.port)]

        command = 'pg_dump {} {}'.format(' '.join(main_args), self.db)
        command += ' > {}'.format(outfile)

        return command

    def send_email(self, outfile):
        _, mail_to = zip(*settings.ADMINS)

        send(
            mail_to,
            settings.DEFAULT_FROM_MAIL,
            template=EmailTemplate.objects.get(name='database_backup_mail'),
            attachments={
                'backup.zip': outfile,
            },
            priority='now',
        )


from .settings import *

DEBUG = False

SECRET_KEY = os.environ.get('DJANGO_SECRET_KEY')

ALLOWED_HOSTS = os.environ.get('ALLOWED_HOSTS')

# Disable template caching
TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
            'loaders' : [
                ('django.template.loaders.cached.Loader',
                                       ['django.template.loaders.filesystem.Loader',
                                        'django.template.loaders.app_directories.Loader'])
            ],
        },
    },
]




# Postgres database
DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': os.environ.get('POSTGRES_DB'),
        'USER': os.environ.get('POSTGRES_USER'),
        'PASSWORD': os.environ.get('POSTGRES_PASSWORD'),
        'HOST': 'database',
        'PORT': '5432',
        'CONN_MAX_AGE': 600,
    }
}

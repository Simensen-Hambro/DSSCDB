from django.apps import AppConfig


class DyeConfig(AppConfig):
    name = 'dye'
    verbose_name = 'Dye data'

    def ready(self):
        from dye import signals

from django.contrib import admin
from .models import Article, Molecule, Spectrum, Performance, Spreadsheet, Contribution, ShortID, AtomicContribution


admin.site.register(Article)
admin.site.register(Molecule)
admin.site.register(Spectrum)
admin.site.register(Performance)
admin.site.register(Spreadsheet)
admin.site.register(Contribution)
admin.site.register(ShortID)
admin.site.register(AtomicContribution)



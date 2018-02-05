from django.contrib import admin
from .models import Article, Molecule, Spectrum, Performance, Spreadsheet, Contribution, AtomicContribution


@admin.register(Molecule)
class EventAdmin(admin.ModelAdmin):
    search_fields = ('inchi',)


admin.site.register(Article)
admin.site.register(Spectrum)
admin.site.register(Performance)
admin.site.register(Spreadsheet)
admin.site.register(Contribution)
admin.site.register(AtomicContribution)

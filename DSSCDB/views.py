from django.shortcuts import render


def index(request):
    context = {}

    return render(request, 'DSSCDB/index.html', context)



from django.shortcuts import render, HttpResponse


def home(request):
    return render(request, "home.html")


def how_it_works(request):
    return render(request, "how-it-works.html")
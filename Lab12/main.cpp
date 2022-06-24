#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const double poczatek = -1.0;
const double koniec = 1.0;
const int N = 11;

double X[N];

double funkcja (double x) {
    return 1.0 / (1.0 + 10.0 * pow(x, 6.0));
}

void wezly_rownoodlegle() {
    for (int i = 0; i < N; ++i) {
        X[i] = poczatek + i * ((koniec - poczatek) / (N - 1.0));
    }
}

void wezly_czebyszewa() {
    double temp;

    for (int i = 0; i < N; ++i) {
        temp = cos(((2.0 * i + 1.0) / (2.0 * N + 2)) * M_PI);
        X[i] = (koniec + poczatek) / 2.0 + ((koniec - poczatek) / 2.0) * temp;
    }
}

void funkcja_dokladna() {
    fstream dokladna;
    dokladna.open("dokladna.txt", fstream::out);

    double i = poczatek;

    while (i <= koniec) {
        dokladna << i << " " << funkcja(i) << "\n";
        i += 0.01;
    }

    dokladna.close();
}

void baza_newtona(string nazwa) {
    fstream plik;
    double C[N][N];

    for (int i = 0; i < N; ++i) {
        C[i][0] = funkcja(X[i]);
    }

    for (int i = 1; i < N; ++i) {
        for (int j = 0; j < N - 1; ++j) {
            C[j][i] = (C[j + 1][i - 1] - C[j][i - 1]) / (X[i + j] - X[j]);
        }
    }

    plik.open(nazwa, fstream::out);

    double a = poczatek;
    double b;

    while (a <= koniec) {
        b = C[0][N - 1];
        for (int i = N - 1; i > 0; --i) {
            b *= (a - X[i - 1]);
            b += C[0][i - 1];
        }
        plik << a << " " << b << "\n";
        a += 0.01;
    }


    plik.close();
}

int main() {
    wezly_rownoodlegle();
    baza_newtona("rownoodlegle.txt");
    wezly_czebyszewa();
    baza_newtona("czebyszewa.txt");
    funkcja_dokladna();
}
#include <iostream>
#include <math.h>
#include <fstream>

const int iteracje = 150;

template <typename T>
T funkcja (T x) {
    return sin(x);
}

template <typename T>
T pochodna_dokladna (T x) {
    return cos(x);
}

template <typename T>
T roznica_progresywna_2pkt (T x, T h) {
    return (funkcja(x + h) - funkcja(x)) / h;
}

template <typename T>
T roznica_wsteczna_2pkt (T x, T h) {
    return (funkcja(x) - funkcja(x - h)) / h;
}

template <typename T>
T roznica_centralna_2pkt (T x, T h) {
    return (funkcja(x + h) - funkcja(x - h)) / (2.0 * h);
}

template <typename T>
T roznica_progresywna_3pkt (T x, T h) {
    return (-3.0 / 2.0 * funkcja(x) + 2.0 * funkcja(x + h) - 1.0 / 2.0 * funkcja(x + h + h)) / h;
}

template <typename T>
T roznica_wsteczna_3pkt (T x, T h) {
    return (1.0 / 2.0 * funkcja(x - h - h) - 2.0 * funkcja(x - h) + 3.0 / 2.0 * funkcja(x)) / h;
}

template <typename T>
void obliczenia (std::string nazwa) {
    T blad[8];
    std::fstream f;
    f.open(nazwa, std::ios::out);

    T h = 0.1; // krok
    T a = 0.0; // poczatek przedzialu
    T b = M_PI_2; // koniec przedzialu
    T c = M_PI_4; //srodek przedzialu

    for (int i = 0; i < iteracje; ++i) {
        blad[0] = log10(h);
        blad[1] = log10(fabs(pochodna_dokladna(a) - roznica_progresywna_2pkt(a, h)));
        blad[2] = log10(fabs(pochodna_dokladna(c) - roznica_progresywna_2pkt(c, h)));
        blad[3] = log10(fabs(pochodna_dokladna(c) - roznica_wsteczna_2pkt(c, h)));
        blad[4] = log10(fabs(pochodna_dokladna(b) - roznica_wsteczna_2pkt(b, h)));
        blad[5] = log10(fabs(pochodna_dokladna(c) - roznica_centralna_2pkt(c, h)));
        blad[6] = log10(fabs(pochodna_dokladna(a) - roznica_progresywna_3pkt(a, h)));
        blad[7] = log10(fabs(pochodna_dokladna(b) - roznica_wsteczna_3pkt(b, h)));

        h /= 1.25;

        for (int j = 0; j < 8; j++) {
            f << blad[j] << " ";
        }
        f << "\n";
    }
}

int main() {
    obliczenia<float>("wynik_float.txt");
    obliczenia<double>("wynik_double.txt");
}
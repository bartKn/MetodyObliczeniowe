#include <iostream>
#include <fstream>
using namespace std;

double taylor(double x) { // x / 1! - x^2 / 2! + x^3 / 3! - x^4 / 4!...
    double suma = x, ostatni = x, znak = -1.0, silnia = 2;

    for (int i = 0; i < 25; ++i) {
        ostatni = ostatni * (x / silnia);
        suma += znak * ostatni;
        silnia += 1.0;
        znak = -znak;
    }
    return suma / x;    // zwracamy / x, zeby miec to zalatwione juz tu
}

double funkcja(double x) {
    double wynik = exp(-x);
    wynik = 1 - wynik;
    wynik /= x;
    return wynik;
}

int main(int argc, char** argv)
{
    cout.precision(20);
    //cout.setf(ios_base::scientific);

    ifstream fin("MO.txt");
    fin.precision(20);

    ofstream fout("wyniki_funkcja.txt");
    fout.precision(20);

    ofstream fout2("wyniki_taylor_25_iteracji.txt");
    fout2.precision(20);

    ofstream fout4("wyniki_z_opisem.txt");
    fout4.precision(20);

    ofstream fout5("wyniki_idelne.txt");
    fout4.precision(20);

    double lg10, x, wynik_dokladny, wynik_obliczony, blad, log_blad;

    double wynik_wzor, blad_wzor, log_blad_wzor;
    double wynik_tylor, blad_taylor, log_blad_taylor;

    cout << "                   lg10(x)|            lg10(blad_wzor)|         lg10(blad_taylor)|\n";

    while (fin >> lg10 >> x >> wynik_dokladny) {
        wynik_wzor = funkcja(x);
        blad_wzor = fabs((wynik_wzor - wynik_dokladny) / wynik_dokladny);
        log_blad_wzor = log10(blad_wzor);
        fout << lg10 << " " << log_blad_wzor << "\n";

        wynik_tylor = taylor(x);
        blad_taylor = fabs((wynik_tylor - wynik_dokladny) / wynik_dokladny);
        log_blad_taylor = log10(blad_taylor);
        fout2 << lg10 << " " << log_blad_taylor << "\n";
        cout.width(25);
        cout << lg10 << " | ";
        cout.width(25);
        cout << log_blad_wzor << " | ";
        cout.width(25);
        cout << log_blad_taylor << "|\n";




        if (log_blad_wzor > log_blad_taylor) {
            fout4 << lg10 << " " << log_blad_taylor << " " << "Taylor\n";
            fout5 << lg10 << " " << log_blad_taylor << "\n";
        } else {
            fout4 << lg10 << " " << log_blad_wzor << " " << "Wzor\n";
            fout5 << lg10 << " " << log_blad_wzor << "\n";
        }

    }

    fin.clear();
    fin.seekg(0, fin.beg);

    ofstream fout3("wyniki_optymalne.txt");
    fout3.precision(20);

    while (fin >> lg10 >> x >> wynik_dokladny) {
        if (lg10 < -0.75) {
            wynik_obliczony = taylor(x);
            blad = fabs((wynik_obliczony - wynik_dokladny) / wynik_dokladny);
            log_blad = log10(blad);
        } else {
            wynik_obliczony = funkcja(x);
            blad = fabs((wynik_obliczony - wynik_dokladny) / wynik_dokladny);
            log_blad = log10(blad);
        }

        fout3 << lg10 << " " << log_blad << "\n";
    }

    fin.close();
    fout.close();
    fout2.close();
    fout3.close();
    fout5.close();

    return 0;
}
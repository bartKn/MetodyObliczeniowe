#include <iostream>
#include <fstream>
#include <cmath>

const double TOLX  = 1e-15;
const double TOLF = 1e-15;
const int MAX_ITER = 30;

using namespace std;

double** utworz_macierz(int n) {
    double** matrix;
    matrix = new double*[n];
    for (int i = 0; i < n; ++i) {
        matrix[i] = new double[n];
    }
    return matrix;
}

void usun_macierz(double** macierz, int n) {
    for (int i = 0; i < n; ++i) {
        delete [] macierz[i];
    }
    delete [] macierz;
}

void wypisz_macierz(double** macierzA, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << macierzA[i][j] << " ";
        }
        cout << endl;
    }
}

void wypisz_wektor(double* wektorB, int n) {
    for (int i = 0; i < n; ++i) {
        cout << wektorB[i] << " ";
    }
    cout << endl;
}

double* reziduum(double** macierzA, double* wektorB, double* wektorX, int n) {
    double suma = 0;
    double* rez = new double[n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            suma += macierzA[i][j] * wektorX[j];
        }
        rez[i] = fabs(suma - wektorB[i]);
        suma = 0;
    }
    return rez;

}

double max(double* wektor, int n) {
    double max = wektor[0];
    for (int i = 1; i < n; ++i) {
        if (wektor[i] > max) {
            max = wektor[i];
        }
    }
    return max;
}

double* estymator(double* wektorXn, double* wektorX, int n) {
    double* est = new double[n];
    for (int i = 0; i < n; ++i) {
        est[i] = fabs(wektorXn[i] - wektorX[i]);
    }
    return est;
}

void metoda_Jacobiego(double** macierzA, double* wektorB, double* wektorX, int n) {

    double suma = 0.0;
    double rez;
    double est;
    double* Xn = new double[n]; // wektor z kolejnymi przybliżeniami

    cout << endl << "1. Metoda Jacobiego:" << endl << endl;
    cout << left << setw(20) <<
    "Nr. iteracji" << setw(20) <<
    "x0" << setw(20) <<
    "x1" << setw(20) <<
    "x2" << setw(20) <<
    "x3" << setw(20) <<
    "Estymator" << setw(20) <<
    "Reziduum" << endl;

    for (int k = 0; k < MAX_ITER; ++k) {

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    suma += macierzA[i][j] * wektorX[j];
                }
            }
            Xn[i] = (wektorB[i] - suma) / macierzA[i][i];
            suma = 0;
        }

        rez = max(reziduum(macierzA, wektorB, Xn, n), n);
        est = max(estymator(Xn, wektorX, n), n);

        cout <<left << setw(20) <<
        k + 1<< setw(20) <<
        Xn[0] << setw(20) <<
        Xn[1] << setw(20) <<
        Xn[2] << setw(20) <<
        Xn[3] << setw(20) <<
        est << setw(20) <<
        rez << endl;

        if (rez <= TOLF && est <= TOLX ) {
            break;
        }

        for (int i = 0; i < n; ++i) {
            wektorX[i] = Xn[i];
        }
    }
}

void metodaGaussaSeidela(double** macierzA, double* wektorB, double* wektorX, int n) {
    double suma = 0.0;
    double rez;
    double est;
    double* Xp = new double[n];

    cout << endl << "2. Metoda Gaussa-Seidel'a:" << endl << endl;
    cout << left << setw(20) <<
         "Nr. iteracji" << setw(20) <<
         "x0" << setw(20) <<
         "x1" << setw(20) <<
         "x2" << setw(20) <<
         "x3" << setw(20) <<
         "Estymator" << setw(20) <<
         "Reziduum" << endl;

    for (int k = 0; k < MAX_ITER; ++k) {

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    suma += macierzA[i][j] * wektorX[j];
                }
            }
            Xp[i] = wektorX[i];
            wektorX[i] = (wektorB[i] - suma) / macierzA[i][i];
            suma = 0;
        }

        rez = max(reziduum(macierzA, wektorB, wektorX, n), n);
        est = max(estymator(wektorX, Xp, n), n);

        cout <<left << setw(20) <<
             k + 1<< setw(20) <<
             wektorX[0] << setw(20) <<
             wektorX[1] << setw(20) <<
             wektorX[2] << setw(20) <<
             wektorX[3] << setw(20) <<
             est << setw(20) <<
             rez << endl;

        if (rez <= TOLF && est <= TOLX ) {
            break;
        }
    }
}

void metoda_SOR(double** macierzA, double* wektorB, double* wektorX, int n, double omega) {
    double suma = 0.0;
    double rez;
    double est;
    double* Xp = new double[n];

    cout << endl << "1. Metoda SOR:" << endl << endl;
    cout << left << setw(20) <<
         "Nr. iteracji" << setw(20) <<
         "x0" << setw(20) <<
         "x1" << setw(20) <<
         "x2" << setw(20) <<
         "x3" << setw(20) <<
         "Estymator" << setw(20) <<
         "Reziduum" << endl;

    for (int k = 0; k < MAX_ITER; ++k) {

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    suma += macierzA[i][j] * wektorX[j];
                }
            }
            Xp[i] = wektorX[i];
            wektorX[i] = (1 - omega) * wektorX[i] + (omega / macierzA[i][i] * (wektorB[i] - suma));
            suma = 0;
        }

        rez = max(reziduum(macierzA, wektorB, wektorX, n), n);
        est = max(estymator(wektorX, Xp, n), n);

        cout <<left << setw(20) <<
             k + 1<< setw(20) <<
             wektorX[0] << setw(20) <<
             wektorX[1] << setw(20) <<
             wektorX[2] << setw(20) <<
             wektorX[3] << setw(20) <<
             est << setw(20) <<
             rez << endl;

        if (rez <= TOLF && est <= TOLX ) {
            break;
        }
    }

    delete[] Xp;
}

int main() {

    ifstream fin("macierz.txt");

    int n = 0;
    fin >> n;

    double** macierzA = utworz_macierz(n);
    double* wektorB = new double[n];
    double* wektorX = new double[n];

    // Wczytanie danych do macierzy A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fin >> macierzA[i][j];
        }
    }

    // Wczytanie danych do wektora B
    // i uzupełnienie wektora X
    for (int i = 0; i < n; ++i) {
        fin >> wektorB[i];
        wektorX[i] = 2;
    }

    metoda_Jacobiego(macierzA, wektorB, wektorX, n);
    for (int i = 0; i < n; ++i) {
        wektorX[i] = 2;
    }
    metodaGaussaSeidela(macierzA, wektorB, wektorX, n);
    for (int i = 0; i < n; ++i) {
        wektorX[i] = 2;
    }
    metoda_SOR(macierzA, wektorB, wektorX, n, 0.5);
    fin.close();

    delete[] wektorX;
    delete[] wektorB;
    usun_macierz(macierzA, n);
    return 0;
}
#include <iostream>
#include <cmath>
#include <fstream>
#include "calerf.h"

using namespace std;

const int D = 1;
const double T_MIN = 0.0;
const double T_MAX = 2.0;
const double A = 6.0 * sqrt(D * T_MAX);
const double X_MIN = 0.0;
const double X_MAX = A;
const double LAMBDA_BEZPOSREDNIE = 0.4;
const double LAMBDA_POSREDNIE = 1.0;
const double H = 0.1;

double rozwiazanie_analityczne (double x, double t) {
    return calerf::ERFC_L(x / (2.0 * sqrt(D * t)));
}

void warunek_poczatkowy (double** macierz, int n, int m) {
    for (int i = 0; i < m; i++) {
        macierz[0][i] = 0.0;
    }
}

void warunki_brzegowe (double** macierz, int n, int m) {
    for (int i = 0; i < n; i++) {
        macierz[i][0] = 1.0;
        macierz[i][m - 1] = 0.0;
    }
}

double** utworz_macierz(int n, int m) {
    double** macierz = new double*[n];

    for (int i = 0; i < n; ++i) {
        macierz[i] = new double[m];
    }
    return macierz;
}

void usun_macierz(double** macierz, int n) {
    for (int i = 0; i < n; ++i) {
        delete [] macierz[i];
    }
    delete [] macierz;
}

void zapisz_macierz(double** macierz, int n, int m, string nazwa_pliku) {
    fstream plik;
    plik.open(nazwa_pliku, fstream::out);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            plik << macierz[i][j] << " ";
        }
        plik << "\n";
    }
    plik.close();
}

void zapisz_wektor(double* wektor, int n, string nazwa_pliku) {
    fstream plik;
    plik.open(nazwa_pliku, fstream::out);

    for (int i = 0; i < n; i++) {
        plik << wektor[i] << "\n";
    }
    plik.close();
}

void zapisz_dwa_wektory(double* wektor1, double* wektor2, int n, string nazwa_pliku) {
    fstream plik;
    plik.open(nazwa_pliku, fstream::out);

    for (int i = 0; i < n; i++) {
        plik << wektor1[i] << " " << wektor2[i] << "\n";
    }
    plik.close();
}

double norma_max (double* wektor, int n) {
    double max = fabs(wektor[0]);

    for (int i = 1; i < n; i++) {
        if (wektor[i] > max) {
            max = wektor[i];
        }
    }
    return max;
}


int wybor_elementu_podstawowego(double** macierz, int n, int j, int* index) {
    int nr_wiersza;

    for (int i = j; i < n - 1; ++i) {
        if (fabs(macierz[index[i]][j]) < fabs(macierz[index[i + 1]][j]))
            nr_wiersza = index[i + 1];
        else
            nr_wiersza = index[i];
    }

    return nr_wiersza;
}

void gauss(double** macierz, double** macierzL, int n, int* index) {

    int nr_wiersza;
    double temp;

    // Przejśćie po diagonali
    for (int i = 0; i < n - 1; ++i) {

        // Sprawdzenie czy 0 na diagonali, jeśli tak konieczny wybór elementu podstawowego
        if (macierz[index[i]][i] == 0.0) {
            nr_wiersza = wybor_elementu_podstawowego(macierz, n, index[i], index);
            // Zamiana wierszy, aktualizujemy dane w tablicy z numerami wierszy
            index[nr_wiersza] = index[i];
            index[i] = nr_wiersza;
        }

        macierzL[index[i]][i] = 1.0;

        // Przejście w dół kolumny, zaczynamy od elementu poniżej diagonali
        for (int j = i + 1; j < n; ++j) {

            // Współczynnik przez który mnożymy przy odejmowaniu
            temp = macierz[index[j]][i] / macierz[index[i]][i];

            // Przejście po wierszu w prawo
            for (int k = i; k < n; ++k) {
                macierz[index[j]][k] = macierz[index[j]][k] - temp * macierz[index[i]][k];
            }
            // W macierzy L umieszczamy współczynniki przez które mnożymy
            macierzL[index[j]][i] = temp;
        }
    }

    // Uzupełnienie ostatniego elementu na diagonali w macierzy L
    macierzL[index[n - 1]][n - 1] = 1.0;
}

void obliczX(double** macierzU, double* wektorX, double* wektorY, int n, int* index) {

    wektorX[index[n - 1]] = wektorY[index[n - 1]] / macierzU[index[n - 1]][n - 1];

    for (int i = n - 1; i >= 0; i--) {
        wektorX[index[i]] = wektorY[index[i]];
        for (int j = i + 1; j < n; j++)
            wektorX[index[i]] = wektorX[index[i]] - macierzU[index[i]][j] * wektorX[index[j]];
        wektorX[index[i]] = wektorX[index[i]] / macierzU[index[i]][i];

    }
}

void obliczY(double** macierzL, double* wektorY, double* wektorB, int n, int* index) {

    wektorY[index[0]] = wektorB[index[0]];

    for (int i = 1; i < n; i++) {
        wektorY[index[i]] = wektorB[index[i]];
        for (int j = 0; j < i; j++) {
            wektorY[index[i]] = wektorY[index[i]] - macierzL[index[i]][j] * wektorY[index[j]];
        }
    }
}

void kmb(double** macierz, int n, int m) {
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m - 1; j++) {
            macierz[i][j] =
                    macierz[i - 1][j] + LAMBDA_BEZPOSREDNIE *
                    (macierz[i - 1][j - 1] - (2.0 * macierz[i - 1][j]) + macierz[i - 1][j + 1]);
        }
    }
}

double** oblicz_kmb(int n, int m) {
    double** macierz_kmb = utworz_macierz(n, m);
    warunek_poczatkowy(macierz_kmb, n, m);
    warunki_brzegowe(macierz_kmb, n, m);
    kmb(macierz_kmb, n, m);
    return macierz_kmb;
}

void thomas_d(double* l, double* d, double* u, int n) {
    for (int i = 1; i < n; i++) {
        d[i] = d[i] - ((l[i - 1] / d[i - 1]) * u[i - 1]);
    }
}

void thomas_bx(double* l, double* d, double* u, double* b, double* x, int n) {
    for (int i = 1; i < n; i++) {
        b[i] = b[i] - ((l[i - 1] / d[i - 1]) * b[i - 1]);
    }

    x[n - 1] = b[n - 1] / d[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
    }
}

void mcn_thomas(double** A, int n, int m) {
    double* l = new double[m];
    double* d = new double[m];
    double* u = new double[m];
    double* b = new double[m];
    double* x = new double[m];

    for (int k = 1; k < n; ++k) {
        d[0] = 1.0;
        u[0] = 0.0;
        b[0] = A[k - 1][0];

        for (int i = 1; i < m - 1; ++i) {
            l[i - 1] = LAMBDA_POSREDNIE / 2.0;
            d[i] = -(1.0 + LAMBDA_POSREDNIE);
            u[i] = LAMBDA_POSREDNIE / 2.0;
            b[i] = -(LAMBDA_POSREDNIE / 2.0 * A[k - 1][i - 1]
                    + (1.0 - LAMBDA_POSREDNIE) * A[k - 1][i]
                    + (LAMBDA_POSREDNIE / 2.0) * A[k - 1][i + 1]);
        }

        l[m - 2] = 0.0;
        d[m - 1] = 1.0;
        b[m - 1] = 0.0;

        thomas_d(l, d, u, m);
        thomas_bx(l, d, u, b, x, m);

        for (int i = 1; i < m - 1; ++i) {
            A[k][i] = x[i];
        }
    }

    delete[] l;
    delete[] d;
    delete[] u;
    delete[] b;
    delete[] x;
}


void mcn_lu(double** A, int n, int m) {
    double* b = new double[m];
    double* y = new double[m];
    double* x = new double[m];
    double* wyniki = new double[m];
    int* index = new int[m];

    double** new_A = utworz_macierz(m, m);
    double** L = utworz_macierz(m, m);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            new_A[i][j] = 0.0;
        }
    }

    for (int i = 0; i < m; ++i) {
        index[i] = i;
        b[i] = 0.0;
        y[i] = 0.0;
        wyniki[i] = 0.0;
    }

    for (int k = 1; k < n; k++) {

        new_A[0][0] = 1.0;
        b[0] = A[k - 1][0];

        for (int i = 1; i < m - 1; i++) {
            new_A[i][i] = -(1.0 + LAMBDA_POSREDNIE);
            new_A[i][i + 1] = LAMBDA_POSREDNIE / 2.0;
            new_A[i][i - 1] = LAMBDA_POSREDNIE / 2.0;
            b[i] = -(LAMBDA_POSREDNIE / 2.0 * A[k - 1][i - 1]
                     + (1.0 - LAMBDA_POSREDNIE) * A[k - 1][i]
                     + (LAMBDA_POSREDNIE / 2.0) * A[k - 1][i + 1]);
        }

        b[m - 1] = 0.0;
        new_A[m - 1][m - 1] = 1.0;

        gauss(new_A, L, m, index);

        obliczY(L, y, b, m, index);

        obliczX(new_A, x, y, m, index);

        for (int i = 0; i < m - 1; ++i) {
            A[k][i] = x[i];
        }
    }
    delete[] b;
    delete[] y;
    delete[] x;
    delete[] wyniki;
    delete[] index;
    usun_macierz(new_A, m);
    usun_macierz(L, m);
}

double** oblicz_mcn_thomas(int n, int m) {
    double** macierz_mcn_thomas = utworz_macierz(n, m);
    warunek_poczatkowy(macierz_mcn_thomas, n, m);
    warunki_brzegowe(macierz_mcn_thomas, n, m);
    mcn_thomas(macierz_mcn_thomas, n, m);
    return macierz_mcn_thomas;
}

double** oblicz_mcn_lu(int n, int m) {
    double** macierz_mcn_lu = utworz_macierz(n, m);
    warunek_poczatkowy(macierz_mcn_lu, n, m);
    warunki_brzegowe(macierz_mcn_lu, n, m);
    mcn_lu(macierz_mcn_lu, n, m);
    return macierz_mcn_lu;
}

double** oblicz_analitycznie(double h, double dt, int n, int m) {
    double** macierz_anal = utworz_macierz(n, m);
    double x = X_MIN;
    double t = T_MIN;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            macierz_anal[i][j] = rozwiazanie_analityczne(x, t);
            x += h;
        }
        x = X_MIN;
        t += dt;
    }
    return macierz_anal;
}

double oblicz_dt(double lambda, double h, double D) {
    return (lambda * h * h) / D;
}

double** oblicz_bledy(double** dokladny, double** obliczony, int n, int m) {
    double** bledy = utworz_macierz(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            bledy[i][j] = fabs(dokladny[i][j] - obliczony[i][j]);
        }
    }
    return bledy;
}

double* oblicz_max_bledy(double** bledy, int n, int m) {
    double* max_bledy = new double[n];
    for (int i = 0; i < n; i++) {
        max_bledy[i] = norma_max(bledy[i], m);
    }
    return max_bledy;
}

double* oblicz_punkty_osi_x(int m) {
    double* punkty = new double[m];
    double x = X_MIN;

    for (int i = 0; i < m; i++) {
        punkty[i] = x;
        x += H;
    }
    return punkty;
}

double* oblicz_punkty_osi_t(double dt, int n) {
    double* punkty = new double[n];
    double t = T_MIN;

    for (int i = 0; i < n; i++) {
        punkty[i] = t;
        t += dt;
    }
    return punkty;
}


void zapis_zad2(double** macierz, double* wektor_kroki, int rozmiar, int pozycja, string nazwa_pliku) {
    double* temp = new double[rozmiar];

    for (int i = 0; i < rozmiar; i++) {
        temp[i] = macierz[pozycja][i];
    }
    zapisz_dwa_wektory(wektor_kroki, temp, rozmiar, nazwa_pliku);
}

void zad2_3_kmb() {
    double dt = oblicz_dt(LAMBDA_BEZPOSREDNIE, H, D);
    int n = (T_MAX - T_MIN) / dt;
    int m = (X_MAX - X_MIN) / H;

    double** rozwiazanie_analityczne;
    double** rozwiazanie_kmb;
    double** bledy;
    double* max_bledy;
    double* punkty_dt;
    double* punkty_x;

    rozwiazanie_analityczne = oblicz_analitycznie(H, dt, n, m);
    zapisz_macierz(rozwiazanie_analityczne, n, m, "rozwiazanie_analityczne.txt");

    rozwiazanie_kmb = oblicz_kmb(n, m);
    zapisz_macierz(rozwiazanie_kmb, n, m, "rozwiazanie_kmb.txt");

    bledy = oblicz_bledy(rozwiazanie_analityczne, rozwiazanie_kmb, n, m);
    max_bledy = oblicz_max_bledy(bledy, n, m);
    zapisz_wektor(max_bledy, n, "max_bledy_kmb.txt");
    zapisz_macierz(bledy, n, m, "bledy_kmb.txt");

    punkty_dt = oblicz_punkty_osi_t(dt, n);
    punkty_x = oblicz_punkty_osi_x(m);

    zapisz_wektor(punkty_dt, n, "odstepy_czasowe_kmb.txt");
    zapisz_wektor(punkty_x, m, "odstepy_x_kmb.txt");

    zapis_zad2(rozwiazanie_kmb, punkty_x, m, 375, "zad2_3_kmb.txt");
    zapis_zad2(rozwiazanie_analityczne, punkty_x, m, 375, "zad2_analitycznie.txt");

    for (int i = 0; i < n; ++i) {
        max_bledy[i] = fabs(max_bledy[i]);
    }

    zapisz_dwa_wektory(punkty_dt, max_bledy, n, "zad3_kmb.txt");

    usun_macierz(rozwiazanie_analityczne, n);
    usun_macierz(rozwiazanie_kmb, n);
    usun_macierz(bledy, n);
    delete[] max_bledy;
    delete[] punkty_dt;
    delete[] punkty_x;
}

void zad2_3_mcn_thomas() {
    double dt = oblicz_dt(LAMBDA_POSREDNIE, H, D);
    int n = (T_MAX - T_MIN) / dt;
    int m = (X_MAX - X_MIN) / H;

    double** rozwiazanie_analityczne;
    double** rozwiazanie_mcn_thomas;
    double** bledy;
    double* max_bledy;
    double* punkty_dt;
    double* punkty_x;

    rozwiazanie_mcn_thomas = oblicz_mcn_thomas(n, m);
    zapisz_macierz(rozwiazanie_mcn_thomas, n, m, "rozwiazanie_mcn_thomas.txt");

    rozwiazanie_analityczne = oblicz_analitycznie(H, dt, n, m);
    zapisz_macierz(rozwiazanie_analityczne, n, m, "rozwiazanie_analityczne.txt");

    bledy = oblicz_bledy(rozwiazanie_analityczne, rozwiazanie_mcn_thomas, n, m);
    max_bledy = oblicz_max_bledy(bledy, n, m);
    zapisz_wektor(max_bledy, n, "max_bledy_mcn_thomas.txt");
    zapisz_macierz(bledy, n, m, "bledy_mcn_thomas.txt");

    punkty_dt = oblicz_punkty_osi_t(dt, n);
    punkty_x = oblicz_punkty_osi_x(m);

    zapisz_wektor(punkty_dt, n, "odstepy_czasowe_mcn_thomas.txt");
    zapisz_wektor(punkty_x, m, "odstepy_x_mcn_thomas.txt");

    zapis_zad2(rozwiazanie_mcn_thomas, punkty_x, m, 150, "zad2_3_mcn_thomas.txt");
    zapis_zad2(rozwiazanie_analityczne, punkty_x, m, 150, "zad2_analitycznie.txt");

    for (int i = 0; i < n; ++i) {
        max_bledy[i] = fabs(max_bledy[i]);
    }

    zapisz_dwa_wektory(punkty_dt, max_bledy, n, "zad3_mcn_t.txt");

    usun_macierz(rozwiazanie_analityczne, n);
    usun_macierz(rozwiazanie_mcn_thomas, n);
    usun_macierz(bledy, n);
    delete[] max_bledy;
    delete[] punkty_dt;
    delete[] punkty_x;
}

void zad2_3_mcn_lu() {
    double dt = oblicz_dt(LAMBDA_POSREDNIE, H, D);
    int n = (T_MAX - T_MIN) / dt;
    int m = (X_MAX - X_MIN) / H;

    double** rozwiazanie_analityczne;
    double** rozwiazanie_mcn_lu;
    double** bledy;
    double* max_bledy;
    double* punkty_dt;
    double* punkty_x;

    rozwiazanie_mcn_lu = oblicz_mcn_lu(n, m);
    zapisz_macierz(rozwiazanie_mcn_lu, n, m, "rozwiazanie_mcn_lu.txt");

    rozwiazanie_analityczne = oblicz_analitycznie(H, dt, n, m);
    zapisz_macierz(rozwiazanie_analityczne, n, m, "rozwiazanie_analityczne.txt");

    bledy = oblicz_bledy(rozwiazanie_analityczne, rozwiazanie_mcn_lu, n, m);
    max_bledy = oblicz_max_bledy(bledy, n, m);
    zapisz_wektor(max_bledy, n, "max_bledy_mcn_lu.txt");
    zapisz_macierz(bledy, n, m, "bledy_mcn_lu.txt");

    punkty_dt = oblicz_punkty_osi_t(dt, n);
    punkty_x = oblicz_punkty_osi_x(m);

    zapisz_wektor(punkty_dt, n, "odstepy_czasowe_mcn_lu.txt");
    zapisz_wektor(punkty_x, m, "odstepy_x_mcn_lu.txt");

    zapis_zad2(rozwiazanie_mcn_lu, punkty_x, m, 150, "zad2_3_mcn_lu.txt");
    zapis_zad2(rozwiazanie_analityczne, punkty_x, m, 150, "zad2_analitycznie.txt");

    for (int i = 0; i < n; ++i) {
        max_bledy[i] = fabs(max_bledy[i]);
    }

    zapisz_dwa_wektory(punkty_dt, max_bledy, n, "zad3_mcn_lu.txt");

    usun_macierz(rozwiazanie_analityczne, n);
    usun_macierz(rozwiazanie_mcn_lu, n);
    usun_macierz(bledy, n);
    delete[] max_bledy;
    delete[] punkty_dt;
    delete[] punkty_x;
}

void zad1() {
    int k = 140;
    double** rozwiazanie_kmb;
    double** rozwiazanie_mcn_t;
    double** rozwiazanie_mcn_lu;
    double** rozwiazanie_analityczne;
    double** bledy;
    double* blad;
    double* wykres_kmb = new double[k];
    double* wykres_mcn_t = new double[k];
    double* wykres_mcn_lu = new double[k];
    double* wykres_h = new double[k];
    double dt;


    double h = 0.25;

    for (int i = 0; i < k; ++i) {
        dt = oblicz_dt(LAMBDA_BEZPOSREDNIE, h, D);
        int n = (T_MAX - T_MIN) / dt;
        int m = (X_MAX - X_MIN) / h;

        rozwiazanie_analityczne = oblicz_analitycznie(h, dt, n, m);
        rozwiazanie_kmb = oblicz_kmb(n, m);
        bledy = oblicz_bledy(rozwiazanie_analityczne, rozwiazanie_kmb, n, m);
        blad = oblicz_max_bledy(bledy, n, m);
        wykres_kmb[i] = log10(fabs(blad[n - 1]));

        dt = oblicz_dt(LAMBDA_POSREDNIE, h, D);
        n = (T_MAX - T_MIN) / dt;
        m = (X_MAX - X_MIN) / h;

        rozwiazanie_analityczne = oblicz_analitycznie(h, dt, n, m);

        rozwiazanie_mcn_t = oblicz_mcn_thomas(n, m);
        bledy = oblicz_bledy(rozwiazanie_analityczne, rozwiazanie_mcn_t, n, m);
        blad = oblicz_max_bledy(bledy, n, m);
        wykres_mcn_t[i] = log10(fabs(blad[n - 1]));

        rozwiazanie_mcn_lu = oblicz_mcn_lu(n, m);
        bledy = oblicz_bledy(rozwiazanie_analityczne, rozwiazanie_mcn_lu, n, m);
        blad = oblicz_max_bledy(bledy, n, m);
        wykres_mcn_lu[i] = log10(fabs(blad[n - 1]));

        wykres_h[i] = log10(h);

        cout << i << endl;
        h /= 1.01;
    }

    zapisz_dwa_wektory(wykres_h, wykres_kmb, k, "zad1_kmb.txt");
    zapisz_dwa_wektory(wykres_h, wykres_mcn_t, k, "zad1_mcn_t.txt");
    zapisz_dwa_wektory(wykres_h, wykres_mcn_lu, k, "zad1_mcn_lu.txt");
}

int main () {
    //zad1();
    //zad2_3_kmb();
    //zad2_3_mcn_thomas();
    //zad2_3_mcn_lu();
}



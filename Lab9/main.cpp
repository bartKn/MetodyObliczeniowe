#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const double P = 1.0, Q = 0.0, R = -4.0;    // Współczynniki zadanego równania

// Z warunków brzegowych
const double alfa = 0.0, beta = 1.0, gama = -1.0;
const double fi = 0.0, psi = 1.0, teta = 0.0;


const double poczatek = 0.0, koniec = 1.0;

double analityczna (double x) {
    return (exp(2.0 - 2.0 * x) - 4.0 * exp(4.0 - x * 2.0) + 4.0 * exp(x * 2.0) - exp(2.0 + 2.0 * x) - x + x * exp(4.0)) / (4.0 - 4 * exp(4.0));
}

void thomas(double *l, double *d, double *u, double *b, double *x, int N) {
    double *r = new double[N];
    double *ni = new double[N];

    ni[0] = d[0];
    r[0] = b[0];

    for (int i = 1; i < N; i++) {
        ni[i] = d[i] - l[i - 1] * (u[i - 1] / ni[i - 1]);
    }

    for (int i = 1; i < N; i++) {
        r[i] = b[i] - l[i - 1] * r[i - 1] / ni[i - 1];
    }

    x[N - 1] = r[N - 1] / ni[N - 1];

    for (int i = N - 2; i >= 0; i--) {
        x[i] = (r[i] - u[i] * x[i + 1]) / ni[i];
    }

    delete[] r;
    delete[] ni;
}


int max_blad (double *blad, int N) {
    double max = blad[0];
    int index = 0;

    for (int i = 0; i < N; i++)
        if (blad[i] > max) {
            max = blad[i];
            index = i;
        }

    return index;
}

double numerow(double h, int N) {
    double *l, *d, *u, *b, *x, *blad;
    double xn = poczatek;

    l = new double[N];
    d = new double[N];
    u = new double[N];
    b = new double[N];
    x = new double[N];
    blad = new double[N];

    u[0] = alfa / h;
    d[0] = beta - alfa / h;
    b[0] = -gama;


    for (int i = 1; i < N - 1; i++) {
        l[i - 1] = P / (h * h) + R / 12.0;
        d[i] = (-2.0 * P) / (h * h) + R * (10.0 / 12.0);
        u[i] = P / (h * h) + R / 12.0;
        b[i] = (xn + i * h - h) / 12.0 + (10.0 / 12.0) * (xn + i * h) + (xn + i * h + h) / 12.0;
    }

    l[N - 2] = -fi / h;
    d[N - 1] = -fi / h + psi;
    b[N - 1] = -teta;


    thomas(l, d, u, b, x, N);

    for (int i = 0; i < N; i++) {
        blad[i] = fabs(x[i] - analityczna(xn));
        xn += h;
    }

    int index = max_blad(blad, N);


    if (N == 1002) {

        fstream numerow;
        numerow.open("numerow.txt", fstream::out);
        numerow << scientific;
        cout.precision(10);


        xn = poczatek;
        for (int i = 0; i < N; i++) {
            numerow << xn << " " << x[i] << "\n";
            xn += h;
        }

        numerow.close();
    }

    delete[] l;
    delete[] d;
    delete[] u;
    delete[] x;
    delete[] b;

    return blad[index];
}

double konwencjonalna(double h, int N) {
    double *l, *d, *u, *b, *x, *blad;
    double xn = poczatek;

    l = new double[N];
    d = new double[N];
    u = new double[N];
    b = new double[N];
    x = new double[N];
    blad = new double[N];


    u[0] = alfa / h;
    d[0] = beta - alfa / h;
    b[0] = -gama;

    for (int i = 1; i < N - 1; i++) {
        l[i - 1] = P / (h * h) - Q / (2.0 * h);
        d[i] = (-2.0 * P) / (h * h) + R;
        u[i] = P / (h * h) + Q / (2.0 * h);
        b[i] = (xn + i * h);
    }

    l[N - 2] = -fi / h;
    d[N - 1] = -fi / h + psi;
    b[N - 1] = -teta;

    thomas(l, d, u, b, x, N);


    for (int i = 0; i < N; i++) {
        blad[i] = fabs(x[i] - analityczna(xn));
        xn += h;
    }

    int index = max_blad(blad, N);

    if (N == 1002) {
        for (int i = 0; i < N; i++) {
            fstream konwencjonalnie, analitycznie;
            konwencjonalnie.open("konwencjonalnie.txt", fstream::out);
            analitycznie.open("analityczna.txt", fstream::out);
            analitycznie << scientific;
            konwencjonalnie << scientific;
            cout.precision(10);


            xn = poczatek;
            for (int i = 0; i < N; i++) {
                konwencjonalnie << xn << " " << x[i] << "\n";
                analitycznie << xn << " " << analityczna(xn) << "\n";
                xn += h;
            }

            analitycznie.close();
            konwencjonalnie.close();
        }
    }

    delete[] l;
    delete[] d;
    delete[] u;
    delete[] x;
    delete[] b;

    return blad[index];
}

int main() {
    double h; // krok
    int N; //ilość iteracji

    fstream bledy;
    bledy.open("bledy.txt", fstream::out);
    bledy << std::scientific;
    cout.precision(10);

    for (N = 2; N < 1000000; N += 50) {
        h = (koniec - poczatek) / (N - 1);
        bledy << log10(h) << " " << log10(numerow(h, N)) << " " << log10(konwencjonalna(h, N)) << "\n";
    }

    bledy.close();

    return 0;
}

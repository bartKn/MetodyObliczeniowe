#include <iostream>

using namespace std;

void wypisz_wektor(double* wektor, int n) {
    for (int i = 0; i < n; ++i) {
        cout << wektor[i] << " ";
    }
    cout << endl;
}

void uzupelnij_wektory(double* U, double* D, double* L, double* B) {
    U[0] = 1.0 / 2.0;       U[1] = 1.0 / 4.0;       U[2] = 1.0 / 6.0;       U[3] = 1.0 / 8.0;       U[4] = 1.0 / 10.0;
    D[0] = 10.0;            D[1] = 20.0;            D[2] = 30.0;            D[3] = 30.0;            D[4] = 20.0;            D[5] = 10.0;
    L[0] = 0.0;             L[1] = 1.0 / 3.0;       L[2] = 1.0 / 5.0;       L[3] = 1.0 / 7.0;       L[4] = 1.0 / 9.0;       L[5] = 1.0 / 11.0;
    B[0] = 31.0;            B[1] = 165.0 / 4.0;     B[2] = 917.0 / 30.0;    B[3] = 851.0 / 28.0;    B[4] = 3637.0 / 90.0;   B[5] = 332.0 / 11.0;
}

void thomas_macierz(double* U, double* D, double* L, double* N, int rozmiar) {
    N[0] = D[0];
    for (int i = 1; i < rozmiar + 1; ++i) {
        N[i] = D[i] - L[i] * U[i - 1] / N[i - 1];
    }
}

void thomas_wektor(double* R, double* B, double* N, double* L, int rozmiar) {
    R[0] = B[0];

    for (int i = 1; i < rozmiar; ++i) {
        R[i] = B[i] - L[i] * R[i - 1] / N[i - 1];
    }
}

void obliczX(double* X, double* N, double* R, double* U, int rozmiar) {
    rozmiar -= 1;
    X[rozmiar] = R[rozmiar] / N[rozmiar];

    for (int i = rozmiar - 1; i >= 0; --i) {
        X[i] = (R[i] - U[i] * X[i + 1]) / N[i];
    }
}

int main() {
    int n = 6;
    double* U = new double[n - 1]; // Przekątna górna
    double* D = new double[n];     // Przekątna główna
    double* L = new double[n - 1]; // Przekątna dolna
    double* B = new double[n];     // Wektor B
    double* X = new double[n];     // Wektor X - wynik
    double* N = new double[n];     // Wektor Eta - nowe wartości na diagonali
    double* R = new double[n];     // Wektor R - nowe wartośći wektora B

    uzupelnij_wektory(U, D, L, B);
    thomas_macierz(U, D, L, N, n);
    thomas_wektor(R, B, N, L, n);
    obliczX(X, N, R, U, n);

    cout << "Wektora Eta:\n";
    wypisz_wektor(N, n);

    cout << "\nWektor R:\n";
    wypisz_wektor(R, n);

    cout << "\nWektor X:\n";
    wypisz_wektor(X, n);



    return 0;
}

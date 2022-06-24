#include <iostream>
#include <cmath>

using namespace std;

typedef double (*funkcja)(double);

double fun1(double x) {             // funkcja_1: sin^2(x / 4) - x
    return sin(x / 4.0) * sin(x / 4.0) - x;
}

double fun1_pochodna(double x) {    // pochodna funkcja_1: 1/4 * sin(x / 2) - 1
    return 0.25 * sin(x / 2.0) - 1.0;
}

double fun1_fi(double x) {          //fi funkcja_1: x = sin^2(x / 4)
    return sin(x / 4.0) * sin(x / 4.0);
}

double fun1_pochodna_fi(double x) { // pochodna fi z funkcja_1: 1/4 * sin(x / 2)
    return 0.25 * sin(x / 2.0);
}

double fun2(double x) {             // funkcja_2: tan(2x) - x - 1
    return tan(2.0 * x) - x - 1.0;
}

double fun2_pochodna(double x) {    // pochodna funkcja_2: (2 / cos^2(2x)) - 1
    return (2.0 / (cos(2.0 * x) * cos(2.0 * x))) - 1.0;
}

double fun2_fi(double x) {          //fi funkcja_2: x = tan(2x) - 1
    return tan(2.0 * x) - 1.0;
}

double fun2_pochodna_fi(double x) { // pochodna fi z funkcja_2: 2 / cos^2(2x)
    return 2.0 / (cos(2.0 * x) * cos(2.0 * x));
}


void picard(funkcja f, funkcja fi, funkcja fi_pochodna, double x, int n_max, double TOLX, double TOLF) {

    cout << "METODA PICARDA\n";

    if (fabs(fi_pochodna(x)) >= 1.0) {          // sprawdzenie zbieżności iteracji, jeśli wartość bezwzględna pochodnej z funkcji fi
        cout << "Brak zbieżności iteracji!\n";  // jest większa lub równa 1, to nie uzyskamy zbieżności
        return;
    }

    double estymator = 0.0;
    double residuum = 0.0;
    double x_nastepny = 0.0;

    cout << "i   Przybliżenie x*\t\t\t\t  Estymator\t\t\t\t\t\tResiduum\n";

    for (int i = 0; i < n_max; ++i) {
        x_nastepny = fi(x);
        estymator = fabs(x_nastepny - x);   // En = Xn - Xn-1
        residuum = fabs(f(x_nastepny));     // f(Xn) - idealnie chcemy uzyskać f(Xn) = 0
        x = x_nastepny;

        cout.width(4); cout << left << i + 1;
        cout.width(30); cout << left << x_nastepny;
        cout.width(30); cout << left << estymator;
        cout.width(30); cout << left << residuum << endl;

        if (estymator <= TOLX && residuum <= TOLF) {    // jeśli spełniony któryś z warunków zakończenia iteracji, przerywamy pętlę
            break;
        }
    }
    cout << "\n\n";
}

void bisekcja(funkcja f, double a, double b, double n_max, double TOLX, double TOLF) {

    cout << "METODA BISEKCJI\n";

    if (f(a) > 0 && f(b) > 0 || f(a) < 0 && f(b) < 0) {         // sprawdzenie założeń początkowych, funkcja musi być różnych znaków
        cout << "W tym przedziale funkcja nie zmienia znaku!\n";  // na końcach przedziału
        return;
    }

    double estymator = 0.0;
    double residuum = 0.0;
    double x = 0.0;         // środek przedziału

    cout << "i\ta\t\t\t\t\t\t\t  b\t\t\t\t\t\t\t\tx\t\t\t\t\t\t\t  Estymator\t\t\t\t\t\tResiduum\n";

    for (int i = 0; i < n_max; ++i) {
        x = (a + b) / 2.0;
        estymator = fabs((b - a) / 2.0);
        residuum = fabs(f(x));

        cout.width(4); cout << left << i + 1;
        cout.width(30); cout << left << a;
        cout.width(30); cout << left << b;
        cout.width(30); cout << left << x;
        cout.width(30); cout << left << estymator;
        cout.width(30); cout << left << residuum << endl;

        if (f(a) > 0 && f(x) < 0 || f(a) < 0 && f(x) > 0) { // wybór nowego przedziału, jeśli f(a) i f(x) są różnych znaków
            b = x;                                          // to b przyjmuje wartość x
        } else {
            a = x;
        }

        if (estymator <= TOLX && residuum <= TOLF) {
            break;
        }
    }
    cout << "\n\n";
}

void newton(funkcja f, funkcja f_pochodna, double x, double n_max, double TOLX, double TOLF) {

    cout << "METODA NEWTONA\n";

    double estymator = 0.0;
    double residuum = 0.0;
    double x_nastepny = 0.0;

    cout << "i   Przybliżenie x*\t\t\t\t  Estymator\t\t\t\t\t\tResiduum\n";

    for (int i = 0; i < n_max; ++i) {
        x_nastepny = x - f(x) / f_pochodna(x);          // Xn+1 = Xn - f(Xn) / f'(Xn)
        estymator = fabs(x_nastepny - x);
        residuum = fabs(f(x_nastepny));
        x = x_nastepny;

        cout.width(4); cout << left << i + 1;
        cout.width(30); cout << left << x_nastepny;
        cout.width(30); cout << left << estymator;
        cout.width(30); cout << left << residuum << endl;

        if (estymator <= TOLX && residuum <= TOLF) {    // jeśli spełniony któryś z warunków zakończenia iteracji, przerywamy pętlę
            break;
        }
    }
    cout << "\n\n";
}

void siecznych(funkcja f, double x0, double x1, int n_max, double TOLX, double TOLF) {

    cout << "METODA SIECZNYCH\n";

    double estymator = 0.0;
    double residuum = 0.0;
    double x_nastepny = 0.0;

    cout << "i   Przybliżenie x*\t\t\t\t  Estymator\t\t\t\t\t\tResiduum\n";

    for (int i = 0; i < n_max; ++i) {
        x_nastepny = x1 - f(x1) / ((f(x1) - f(x0)) / (x1 - x0)); // Xn+2 = Xn+1 - f(Xn+1) / [(f(Xn+1) - f(Xn)) / (Xn+1 - Xn)]
        estymator = fabs(x_nastepny - x1);
        residuum = fabs(f(x_nastepny));
        x0 = x1;
        x1 = x_nastepny;

        cout.width(4); cout << left << i + 1;
        cout.width(30); cout << left << x_nastepny;
        cout.width(30); cout << left << estymator;
        cout.width(30); cout << left << residuum << endl;

        if (estymator <= TOLX && residuum <= TOLF) {    // jeśli spełniony któryś z warunków zakończenia iteracji, przerywamy pętlę
            break;
        }
    }
    cout << "\n\n";
}



int main() {
    //cout.setf(ios::scientific);
    cout.precision(20);

    double TOL = 1e-15;
    int MAX_N = 100;

    cout << "Równanie 1. sin^2(x / 4) - x = 0\n"; // x = 0
    picard(fun1, fun1_fi, fun1_pochodna_fi, 1, MAX_N, TOL, TOL);
    bisekcja(fun1, -1.5, 1.9, MAX_N, TOL, TOL);
    newton(fun1, fun1_pochodna, -1.2, MAX_N, TOL, TOL);
    siecznych(fun1, -1.5, 1.6, MAX_N, TOL, TOL);

    cout << "Równanie 2. tan(x / 2) - x - 1 = 0\n"; // x = -1.951, x = 0.49, x = 2.205
    picard(fun2, fun2_fi, fun2_pochodna_fi, 0.49, MAX_N, TOL, TOL);
    bisekcja(fun2, 0.3, 0.6, MAX_N, TOL, TOL);
    newton(fun2, fun2_pochodna, 0.3, MAX_N, TOL, TOL);
    siecznych(fun2, 0.3, 0.6, MAX_N, TOL, TOL);
}


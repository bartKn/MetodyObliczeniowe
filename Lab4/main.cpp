#include <iostream>
#include <math.h>

using namespace std;

double fun0(double x, double y, double z) { // x^2 + y^2 + z^2 - 2
    return x * x + y * y + z * z - 2.0;
}

double fun1(double x, double y) { // x^2 + y^2 -1
    return x * x + y * y - 1.0;
}

double fun2(double x, double y) { // x^2 - y
    return x * x - y;
}

double delta0(double x, double y, double z) { // [ x^2 - y^2 + 2 * x^2 * y - 1 ] / [ 2*x * (1 + 2*y) ]
    return (x * x - y * y + 2.0 * x * x * y - 1.0) / (2.0 * x * (1.0 + 2.0 * y));
}

double delta1(double x, double y, double z) { // [ y^2 + y -1 ] / [2*y + 1]
    return (y * y + y - 1.0) / (2.0 * y + 1.0);
}

double delta2(double x, double y, double z) { // [ z^2 - 1 ] / [ 2*z ]
    return (z * z - 1.0) / (2 * z);
}


double max(double a, double b, double c) {
    double max = a;
    if (b > max) max = b;
    if (c > max) max = c;
    return max;
}

void uogolnionaMetodaNewtona(double x, double y, double z, int n_max, double TOLX, double TOLF) {
    // wejsciowe x, y, z to elementy wektora niewiadomych
    double x_nastepny= 0.0;
    double y_nastepny = 0.0;
    double z_nastepny = 0.0;
    double estymator = 0.0;
    double residuum = 0.0;
    double delta[3];        // wektor delta = J^1 * f, odwrotność jakobianu * wektor z wartościami funkcji
    double f[3];            // wektor na wartości funkcji w kolejnym przybliżeniu

    cout << "i\tx\t\t\t\t\t\t\t  y\t\t\t\t\t\t\t\tz\t\t\t\t\t\t\t  Estymator\t\t\t\t\t\tResiduum\n";

    for (int i = 0; i < n_max; ++i) {

        delta[0] = delta0(x, y, z);
        delta[1] = delta1(x, y, z);
        delta[2] = delta2(x, y, z);

        // wyliczenie następnych przybliżeń, Xn+1 = Xn - delta
        x_nastepny = x - delta[0];
        y_nastepny = y - delta[1];
        z_nastepny = z - delta[2];

        //obliczenie wartości funkcji dla kolejnych przybliżeń, na ich podstawie obliczamy residuum
        f[0] = fun0(x_nastepny, y_nastepny, z_nastepny);
        f[1] = fun1(x_nastepny, y_nastepny);
        f[2] = fun2(x_nastepny, y_nastepny);

        estymator = max(fabs(x_nastepny - x), fabs(y_nastepny - y), fabs(z_nastepny - z));
        residuum = max(fabs(f[0]), fabs(f[1]), fabs(f[2]));



        x = x_nastepny;
        y = y_nastepny;
        z = z_nastepny;

        cout.width(4); cout << left << i + 1;
        cout.width(30); cout << left << x;
        cout.width(30); cout << left << y;
        cout.width(30); cout << left << z;
        cout.width(30); cout << left << estymator;
        cout.width(30); cout << left << residuum << endl;

        if (estymator <= TOLX && residuum <= TOLF) {
            break;
        }
    }
}

int main() {
    cout.precision(20);
    int N_MAX = 50;
    double TOL = 1e-15;
    uogolnionaMetodaNewtona(-1, 3, -2, N_MAX, TOL, TOL);

    return 0;
}


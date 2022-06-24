#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double analityczne(double t) {
    return 1 - exp(-10.0 * (t + atan(t)));
}

double BME(double t, double delta_t) {
    double yk = 0; // z warunku początkowego
    double tk = 0;

    while (tk < t) {

        yk = yk - ((10.0 * tk * tk + 20.0) / (tk * tk + 1.0)) * (yk - 1.0) * delta_t;

        tk += delta_t;
    }

    return yk;
}

double BME_blad(double delta_t) {
    double max_blad = 0.0;
    double blad;
    double tk = 0.0;
    double yk = 0.0;

    while (tk < 1) {
        blad = fabs(analityczne(tk) - yk);

        if (blad > max_blad) {
            max_blad = blad;
        }

        yk = yk - ((10.0 * tk * tk + 20.0) / (tk * tk + 1.0)) * (yk - 1.0) * delta_t;

        tk += delta_t;
    }

    return max_blad;
}

double PME(double t, double delta_t) {
    double yk = 0.0;
    double tk = 0.0;
    double temp;

    while (tk < t) {

        temp = ((10.0 * (tk + delta_t) * (tk + delta_t) + 20.0) / ((tk + delta_t) * (tk + delta_t) + 1.0)) * delta_t;
        yk = (yk + temp) / (1 + temp);

        tk += delta_t;
    }

    return yk;
}

double PME_blad(double delta_t) {
    double max_blad = 0.0;
    double blad;
    double tk = 0.0;
    double yk = 0.0;
    double temp;

    while (tk < 1.0) {
        blad = fabs(analityczne(tk) - yk);

        if (blad > max_blad) {
            max_blad = blad;
        }

        temp = ((10.0 * (tk + delta_t) * (tk + delta_t) + 20.0) / ((tk + delta_t) * (tk + delta_t) + 1.0)) * delta_t;
        yk = (yk + temp) / (1 + temp);

        tk += delta_t;
    }
    return max_blad;
}

double PMT(double t, double delta_t) {
    double yk = 0.0;
    double tk = 0.0;
    double temp1;
    double temp2;

    while (tk < t) {

        temp1 = (10.0 * tk * tk + 20.0) / (tk * tk + 1.0);
        temp2 = (10.0 * (tk + delta_t) * (tk + delta_t) + 20.0) / ((tk + delta_t) * (tk + delta_t) + 1.0);
        yk = ((-delta_t / 2.0) * (temp1 * (yk - 1.0) - temp2) + yk) / (1.0 + (delta_t / 2.0) * temp2);

        tk += delta_t;
    }

    return yk;
}

double PMT_blad(double delta_t) {
    double max_blad = 0.0;
    double blad;
    double tk = 0.0;
    double yk = 0.0;
    double temp1;
    double temp2;

    while (tk < 1.0) {
        blad = fabs(analityczne(tk) - yk);

        if (blad > max_blad) {
            max_blad = blad;
        }

        temp1 = (10.0 * tk * tk + 20.0) / (tk * tk + 1.0);
        temp2 = (10.0 * (tk + delta_t) * (tk + delta_t) + 20.0) / ((tk + delta_t) * (tk + delta_t) + 1.0);
        yk = ((-delta_t / 2.0) * (temp1 * (yk - 1.0) - temp2) + yk) / (1.0 + (delta_t / 2.0) * temp2);

        tk += delta_t;
    }
    return max_blad;
}


int main() {
    fstream wynikAnalitycznie, wynikiBME, wynikiPME, wynikiPMT, wynikiBlad;

    wynikAnalitycznie.open("wynikAnalitycznie.txt", fstream::out);
    wynikiBME.open("wynikiBME.txt", fstream::out);
    wynikiPME.open("wynikiPME.txt", fstream::out);
    wynikiPMT.open("wynikiPMT.txt", fstream::out);
    wynikiBlad.open("wynikiBlad.txt", fstream::out);

    wynikAnalitycznie << scientific;
    wynikiBME << scientific;
    wynikiPME << scientific;
    wynikiPMT << scientific;
    wynikiBlad << scientific;

    cout.precision(10);

    double delta_t = 0.01;

    double t = 0.0;

    while (t <= 2) {
        wynikAnalitycznie << t << " " << analityczne(t) << "\n";
        wynikiBME << t << " " << BME(t, delta_t) << "\n";
        wynikiPME << t << " " << PME(t, delta_t) << "\n";
        wynikiPMT << t << " " << PMT(t, delta_t) << "\n";

        t += delta_t;
    }

    delta_t = 0.1;

    while (delta_t > 1e-10) {

        wynikiBlad << log10(delta_t) << " " << log10(BME_blad(delta_t)) << " " << log10(PME_blad(delta_t)) << " " << log10(PMT_blad(delta_t)) << "\n";

        delta_t /= 1.5;
    }

    wynikAnalitycznie.close();
    wynikiBME.close();
    wynikiPME.close();
    wynikiPMT.close();
    wynikiBlad.close();

}
#include <iostream>

using namespace std;

/*
 * Epsilon maszynowy - najmiejsza liczba która dodana do jedności daje wartość większą od 1
 *
 * Precyzja arytmetyki to 2 ^ -(t+1)
 *
 * Epsilon jest równy podwojonej precyzji arytmetyki
 *
 *
 * t dla float = 23
 *
 * t dla double = 52
 */

int main() {

    float eps_float = 1.0f;
    float next_eps_float;
    int t = 0;

    do {
        eps_float /= 2.0f;
        t += 1;
        next_eps_float = eps_float / 2.0f;
        next_eps_float += 1.0f;
    } while (next_eps_float > 1.0f);


    /*while( ((eps_float / 2.0f) + 1.0f) > 1.0f ) {
        eps_float /= 2.0f;
        t += 1;
    }*/


    cout << "Epsilon maszynowy dla float = " << eps_float << endl;
    cout << "Liczba bitow mantysy dla float = " << t << endl << endl;


    double eps_double = 1.0;
    double next_eps_double;
    t = 0;

    do {
        eps_double /= 2.0;
        t += 1;
        next_eps_double = eps_double / 2.0;
        next_eps_double += 1.0;
    } while (next_eps_double > 1.0);

   /* while( ((eps_double / 2.0) + 1.0) > 1.0 ) {
        eps_double /= 2.0;
        t += 1;
    }*/

    cout << "Epsilon maszynowy dla double = " << eps_double << endl;
    cout << "Liczba bitow mantysy dla double = " << t << endl;

    return 0;
}


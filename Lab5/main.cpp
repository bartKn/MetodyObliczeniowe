#include <iostream>
#include <fstream>
#include <cmath>

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

void wypisz_macierz(double** macierzA, int n, int* index) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << macierzA[index[i]][j] << " ";
        }
        cout << endl;
    }
}

void wypisz_wektor(double* wektorB, int n, int* index) {
    for (int i = 0; i < n; ++i) {
        cout << wektorB[index[i]] << " ";
    }
    cout << endl;
}

// Wybór elementu podstawowego, szukamy wartości największej w danej kolumnie poniżej komórki [j][j],
// zwraca numer wiersza w którym znajduje się wartosć największa
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


int main() {
    ifstream fin("macierz.txt");

    int n = 0;
    fin >> n;

    double** macierzA = utworz_macierz(n);
    double** macierzL = utworz_macierz(n);
    double* wektorB = new double[n];
    double* wektorX = new double[n];
    double* wektorY = new double[n];

    // Wczytanie danych do macierzy A
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fin >> macierzA[i][j];
        }
    }

    // Wczytanie danych do wektora B
    for (int i = 0; i < n; ++i) {
        fin >> wektorB[i];
    }

    // Utworzenie tablicy z indeksami
    int* indeks;
    indeks = new int[n];
    for (int i = 0; i < n; ++i) {
        indeks[i] = i;
    }

    cout << "Macierz A:\n";
    wypisz_macierz(macierzA, n, indeks);
    cout << "\nWektor B:\n";
    wypisz_wektor(wektorB, n, indeks);

    gauss(macierzA, macierzL, n, indeks);

    cout << "\nMacierz L:\n";
    wypisz_macierz(macierzL, n, indeks);

    cout << "\nMacierz U:\n";
    wypisz_macierz(macierzA, n, indeks);

    obliczY(macierzL, wektorY, wektorB, n, indeks);
    cout << "\nWektor Y:\n";
    wypisz_wektor(wektorY, n, indeks);

    obliczX(macierzA, wektorX, wektorY, n, indeks);
    cout << "\nWektor X:\n";
    wypisz_wektor(wektorX, n, indeks);

    ofstream fout("wynik.txt");

    fout << n << "\n\n";

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << macierzL[indeks[i]][j] << " ";
        }
        fout << "\n";
    }
    fout << "\n";

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << macierzA[indeks[i]][j] << " ";
        }
        fout << "\n";
    }
    fout << "\n";

    for (int i = 0; i < n; ++i) {
        fout << wektorX[indeks[i]] << " ";
    }



    delete[] wektorB;
    delete[] wektorX;
    delete[] wektorY;
    usun_macierz(macierzA, n);
    usun_macierz(macierzL, n);
    fin.close();
    fout.close();

    return 0;
}

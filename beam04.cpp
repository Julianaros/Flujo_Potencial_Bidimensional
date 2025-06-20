#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>  // para std::max

// Parámetros del dominio
const int Nxmax = 70;
const int Nymax = 20;
const int IL = 10;
const int H = 8;
const int T = 8;
const double h = 1.0;
const double V0 = 1.0;
const double omega = 0.1;
const double tol = 1e-5;

double u[Nxmax + 1][Nymax + 1] = {0};  // Función de corriente (psi)

// Condiciones de frontera
void borders() {
    for (int i = 0; i <= Nxmax; ++i)
        for (int j = 0; j <= Nymax; ++j)
            u[i][j] = j * V0;

    for (int i = 0; i <= Nxmax; ++i)
        u[i][Nymax] = Nymax * V0;  // corregido: no extrapolar fuera de rango

    for (int i = 0; i <= Nxmax; ++i)
        if (i <= IL || i >= IL + T)
            u[i][0] = 0.0;

    for (int j = 0; j <= Nymax; ++j)
        u[0][j] = j * V0;

    for (int j = 1; j < Nymax; ++j)
        u[Nxmax][j] = u[Nxmax - 1][j];
}

// Zona bloqueada (la viga)
void block_beam_zone() {
    for (int i = IL; i <= IL + T && i <= Nxmax; ++i)
        for (int j = 1; j <= H && j <= Nymax; ++j)  // j=1 evita sobreescribir frontera j=0
            u[i][j] = 0.0;
}

// Relajación con criterio de convergencia
void relax_until_converge() {
    int iter = 0;
    double max_diff;

    do {
        max_diff = 0.0;
        block_beam_zone();

        for (int i = 1; i < Nxmax; ++i) {
            for (int j = 1; j < Nymax; ++j) {
                if (i >= IL && i <= IL + T && j <= H)
                    continue;

                double old = u[i][j];
                double r1 = omega * ((u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) / 4.0 - old);
                u[i][j] += r1;
                max_diff = std::max(max_diff, std::abs(r1));
            }
        }

        ++iter;

    } while (max_diff > tol);

    std::cout << "Convergencia alcanzada en " << iter << " iteraciones. "
              << "Error máximo final: " << max_diff << "\n";
}


// Normalización
void normalize() {
    for (int i = 0; i <= Nxmax; ++i)
        for (int j = 0; j <= Nymax; ++j)
            u[i][j] /= (V0 * h);
}

// Exportar función de corriente
void export_streamfunction(const std::string& filename) {
    std::ofstream file(filename);
    for (int i = 0; i < Nxmax; ++i) {
        for (int j = 0; j < Nymax; ++j)
            file << i << " " << j << " " << u[i][j] << "\n";
        file << "\n";
    }
    file.close();
}

// Exportar campo de velocidades
void export_velocity_field(const std::string& filename) {
    double max_mag = 0.0;

    for (int i = 2; i < Nxmax - 2; ++i) {
        for (int j = 2; j < Nymax - 2; ++j) {
            if (i >= IL && i <= IL + T && j <= H)
                continue;

            double dudY = (u[i][j+1] - u[i][j-1]) / (2 * h);
            double dvdX = -(u[i+1][j] - u[i-1][j]) / (2 * h);
            double mag = std::sqrt(dudY * dudY + dvdX * dvdX);
            if (mag > max_mag)
                max_mag = mag;
        }
    }

    std::ofstream file(filename);
    for (int i = 2; i < Nxmax - 2; ++i) {
        for (int j = 2; j < Nymax - 2; ++j) {
            if (i >= IL && i <= IL + T && j <= H)
                continue;

            double dudY = (u[i][j+1] - u[i][j-1]) / (2 * h);
            double dvdX = -(u[i+1][j] - u[i-1][j]) / (2 * h);

            if (max_mag > 1e-8) {
                dudY /= max_mag;
                dvdX /= max_mag;
            }

            file << i << " " << j << " " << dudY << " " << dvdX << "\n";
        }
        file << "\n";
    }

    file.close();
}

// Función principal
int main() {
    std::cout << "Calculando flujo con relajación hasta convergencia...\n";
    borders();
    relax_until_converge();
    normalize();
    export_streamfunction("streamfunction.dat");
    export_velocity_field("velocity_field.dat");
    std::cout << "Exportación completa: streamfunction.dat & velocity_field.dat\n";
    return 0;
}

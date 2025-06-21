#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

// Parámetros del dominio
const int Nxmax = 70;
const int Nymax = 20;
const int IL = 10;      // Inicio de la viga
const int H = 8;        // Altura de la viga
const int T = 8;        // Longitud de la viga
const double h = 1.0;   // Espaciado de la malla
const double V0 = 1.0;  // Velocidad de entrada
const double omega = 0.1; // Parámetro de sobre-relajación
const double nu = 1.0;  // Viscosidad cinemática
const double tol = 1e-5; // Tolerancia de convergencia


// Variables principales
double u[Nxmax + 1][Nymax + 1] = {0};  // Función de corriente (stream function)
double w[Nxmax + 1][Nymax + 1] = {0};  // Vorticidad
double R;  // Número de Reynolds de malla

// Inicialización de condiciones de frontera
void borders() {
    // Inicializar función de corriente y vorticidad
    for (int i = 0; i <= Nxmax; ++i) {
        for (int j = 0; j <= Nymax; ++j) {
            w[i][j] = 0.0;          // Inicializar vorticidad
            u[i][j] = j * V0;       // Flujo libre inicial
        }
    }
    
    // Superficie del fluido (condición de flujo libre)
    for (int i = 0; i <= Nxmax; ++i) {
        u[i][Nymax] = u[i][Nymax-1] + V0 * h;
        if (Nymax > 1) w[i][Nymax-1] = 0.0;
    }
    
    // Entrada (inlet) - flujo horizontal uniforme
    for (int j = 0; j <= Nymax; ++j) {
        u[1][j] = u[0][j];  // du/dx = 0
        w[0][j] = 0.0;      // sin rotación
    }
    
    // Línea central y viga
    for (int i = 0; i <= Nxmax; ++i) {
        if (i <= IL || i >= IL + T) {
            u[i][0] = 0.0;      // Línea central: u = 0
            w[i][0] = 0.0;      // Sin vorticidad en línea central
        }
    }
    
    // Salida (outlet) - condiciones de gradiente cero
    for (int j = 1; j < Nymax; ++j) {
        w[Nxmax][j] = w[Nxmax-1][j];  // dw/dx = 0
        u[Nxmax][j] = u[Nxmax-1][j];  // du/dx = 0
    }
}

// Condiciones de frontera específicas para la viga
void beam_boundaries() {
    // Lados de la viga (frente y atrás)
    for (int j = 0; j <= H && j <= Nymax; ++j) {
        if (IL > 0) {
            w[IL][j] = -2.0 * u[IL-1][j] / (h * h);     // Frente de la viga
        }
        if (IL + T < Nxmax) {
            w[IL+T][j] = -2.0 * u[IL+T+1][j] / (h * h); // Atrás de la viga
        }
    }
    
    // Parte superior de la viga
    for (int i = IL; i <= IL + T && i <= Nxmax; ++i) {
        if (H < Nymax) {
            w[i][H-1] = -2.0 * u[i][H] / (h * h);
        }
    }
    
    // Función de corriente en la viga (u = 0)
    for (int i = IL; i <= IL + T && i <= Nxmax; ++i) {
        for (int j = 0; j <= H && j <= Nymax; ++j) {
            u[IL][j] = 0.0;     // Frente
            u[IL+T][j] = 0.0;   // Atrás
            u[i][H] = 0.0;      // Superior
        }
    }
}

// Función de relajación para función de corriente
void relax_stream_function() {
    for (int i = 1; i < Nxmax; ++i) {
        for (int j = 1; j < Nymax; ++j) {
            // Saltear puntos dentro de la viga
            if (i >= IL && i <= IL + T && j <= H) continue;
            
            double old_u = u[i][j];
            double new_u = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] + h*h*w[i][j]);
            double r1 = omega * (new_u - old_u);
            u[i][j] += r1;
        }
    }
}

// Función de relajación para vorticidad
void relax_vorticity() {
    for (int i = 1; i < Nxmax; ++i) {
        for (int j = 1; j < Nymax; ++j) {
            // Saltear puntos dentro de la viga
            if (i >= IL && i <= IL + T && j <= H) continue;
            
            double old_w = w[i][j];
            
            // Términos de la ecuación de vorticidad
            double a1 = w[i+1][j] + w[i-1][j] + w[i][j+1] + w[i][j-1];
            double a2 = (u[i][j+1] - u[i][j-1]) * (w[i+1][j] - w[i-1][j]);
            double a3 = (u[i+1][j] - u[i-1][j]) * (w[i][j+1] - w[i][j-1]);
            
            double new_w = 0.25 * (a1 - (R/4.0) * (a2 - a3));
            double r2 = omega * (new_w - old_w);
            w[i][j] += r2;
        }
    }
}

// Función principal de relajación hasta convergencia
void relax_until_converge() {
    int iter = 0;
    double max_diff_u, max_diff_w;
    
    // Calcular número de Reynolds de malla
    R = V0 * h / nu;
    std::cout << "Número de Reynolds de malla R = " << R << std::endl;
    
    do {
        max_diff_u = 0.0;
        max_diff_w = 0.0;
        
        // Aplicar condiciones de frontera
        beam_boundaries();
        
        // Guardar valores anteriores para calcular diferencias
        double u_old[Nxmax + 1][Nymax + 1];
        double w_old[Nxmax + 1][Nymax + 1];
        
        for (int i = 0; i <= Nxmax; ++i) {
            for (int j = 0; j <= Nymax; ++j) {
                u_old[i][j] = u[i][j];
                w_old[i][j] = w[i][j];
            }
        }
        
        // Relajar función de corriente
        relax_stream_function();
        
        // Relajar vorticidad
        relax_vorticity();
        
        // Calcular diferencias máximas
        for (int i = 1; i < Nxmax; ++i) {
            for (int j = 1; j < Nymax; ++j) {
                if (i >= IL && i <= IL + T && j <= H) continue;
                
                max_diff_u = std::max(max_diff_u, std::abs(u[i][j] - u_old[i][j]));
                max_diff_w = std::max(max_diff_w, std::abs(w[i][j] - w_old[i][j]));
            }
        }
        
        ++iter;
        
        
    } while ((max_diff_u > tol || max_diff_w > tol) && iter < 100000);
    
    std::cout << "Convergencia alcanzada en " << iter << " iteraciones." << std::endl;
    std::cout << "Error final: u = " << max_diff_u << ", w = " << max_diff_w << std::endl;
}

// Normalización
void normalize() {
    for (int i = 0; i <= Nxmax; ++i) {
        for (int j = 0; j <= Nymax; ++j) {
            u[i][j] /= (V0 * h);
        }
    }
}

// Exportar función de corriente
void export_streamfunction(const std::string& filename) {
    std::ofstream file(filename);
    file << "# Función de corriente - Formato: x y psi" << std::endl;
    file << "# Parámetros: Nx=" << Nxmax << " Ny=" << Nymax << " h=" << h << std::endl;
    file << std::fixed << std::setprecision(6);
    
    for (int i = 0; i < Nxmax; ++i) {
        for (int j = 0; j < Nymax; ++j) {
            // Usar coordenadas físicas en lugar de índices
            double x = i * h;
            double y = j * h;
            file << x << " " << y << " " << u[i][j] << std::endl;
        }
    }
    file.close();
    std::cout << "Función de corriente exportada a: " << filename << std::endl;
}

// Exportar vorticidad
void export_vorticity(const std::string& filename) {
    std::ofstream file(filename);
    file << "# Vorticidad - Formato: x y omega" << std::endl;
    file << "# Parámetros: Nx=" << Nxmax << " Ny=" << Nymax << " h=" << h << std::endl;
    file << std::fixed << std::setprecision(6);
    
    for (int i = 0; i < Nxmax; ++i) {
        for (int j = 0; j < Nymax; ++j) {
            double x = i * h;
            double y = j * h;
            file << x << " " << y << " " << w[i][j] << std::endl;
        }
    }
    file.close();
    std::cout << "Vorticidad exportada a: " << filename << std::endl;
}

// Exportar campo de velocidades
void export_velocity_field(const std::string& filename) {
    std::ofstream file(filename);
    file << "# Campo de velocidades - Formato: x y vx vy velocidad_magnitud" << std::endl;
    file << "# Parámetros: Nx=" << Nxmax << " Ny=" << Nymax << " h=" << h << std::endl;
    file << std::fixed << std::setprecision(6);
    
    for (int i = 1; i < Nxmax - 1; ++i) {
        for (int j = 1; j < Nymax - 1; ++j) {
            // Saltear puntos dentro de la viga
            if (i >= IL && i <= IL + T && j <= H) continue;
            
            double x = i * h;
            double y = j * h;
            
            // Calcular componentes de velocidad: vx = du/dy, vy = -du/dx
            double vx = (u[i][j+1] - u[i][j-1]) / (2.0 * h);
            double vy = -(u[i+1][j] - u[i-1][j]) / (2.0 * h);
            double v_mag = std::sqrt(vx*vx + vy*vy);
            
            file << x << " " << y << " " << vx << " " << vy << " " << v_mag << std::endl;
        }
    }
    file.close();
    std::cout << "Campo de velocidades exportado a: " << filename << std::endl;
}

// Función principal
int main() {
    std::cout << "=== Solver de Navier-Stokes (Forma de Vorticidad) ===" << std::endl;
    std::cout << "Parámetros:" << std::endl;
    std::cout << "  Malla: " << Nxmax << " x " << Nymax << std::endl;
    std::cout << "  Viga: posición x=[" << IL << "," << IL+T << "], altura=" << H << std::endl;
    std::cout << "  V0 = " << V0 << ", nu = " << nu << ", omega = " << omega << std::endl;
    std::cout << std::endl;
    
    // Inicializar condiciones de frontera
    borders();
    
    // Resolver hasta convergencia
    relax_until_converge();
    
    // Normalizar resultados
    normalize();
    
    // Exportar resultados
    export_streamfunction("streamfunction.dat");
    export_vorticity("vorticity.dat");
    export_velocity_field("velocity_field.dat");
    
    std::cout << std::endl << "¡Simulación completada exitosamente!" << std::endl;
    std::cout << "Archivos generados:" << std::endl;
    std::cout << "  - streamfunction.dat: Función de corriente" << std::endl;
    std::cout << "  - vorticity.dat: Campo de vorticidad" << std::endl;
    std::cout << "  - velocity_field.dat: Campo de velocidades" << std::endl;
    
    return 0;
}
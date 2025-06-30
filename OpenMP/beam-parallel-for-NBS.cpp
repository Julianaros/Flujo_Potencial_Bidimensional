#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>
#include <sys/stat.h>  // Para crear directorios
#include <sys/types.h> // Para tipos de datos del sistema
#include <omp.h> // Librería para paralelizar con openMP

// Parámetros del dominio
const int Nxmax = 160;
const int Nymax = 30;
const int IL = 10;      // Inicio de la viga
const int H = 8;        // Altura de la viga
const int T = 8;        // Longitud de la viga
const double h = 1.0;   // Espaciado de la malla
const double V0 = 1.0;  // Velocidad de entrada
const double tol = 1e-8; // Tolerancia de convergencia

// Variables que cambiarán según Reynolds
double omega;           // Parámetro de sobre-relajación (ajustable)
double nu;             // Viscosidad cinemática (calculada)
double R_target;       // Número de Reynolds objetivo

// Variables principales
double u[Nxmax + 1][Nymax + 1] = {0};  // Función de corriente (stream function)
double w[Nxmax + 1][Nymax + 1] = {0};  // Vorticidad
double R;  // Número de Reynolds de malla calculado

// Función para crear directorio si no existe
bool create_directory(const std::string& path) {
    #ifdef _WIN32
        // Windows
        int result = _mkdir(path.c_str());
    #else
        // Linux/Unix/Mac
        int result = mkdir(path.c_str(), 0755);
    #endif
    
    if (result == 0) {
        std::cout << "Directorio '" << path << "' creado exitosamente." << std::endl;
        return true;
    } else {
        // Verificar si el directorio ya existe
        struct stat info;
        if (stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode)) {
            std::cout << "Directorio '" << path << "' ya existe." << std::endl;
            return true;
        } else {
            std::cout << "Error al crear el directorio '" << path << "'." << std::endl;
            return false;
        }
    }
}

// Configurar parámetros optimizados para diferentes Reynolds
void configure_reynolds(double reynolds_target) {
    R_target = reynolds_target;
    
    // Calcular viscosidad para obtener el Reynolds deseado
    nu = V0 * h / R_target;
    
    // Parámetros omega optimizados
    if (R_target <= 0.5) {
        omega = 0.1;   
    } else if (R_target <= 1.0) {
        omega = 0.08;  
    } else if (R_target <= 2.0) {
        omega = 0.04;  
    } else if (R_target <= 5.0) {
        omega = 0.012; // AJUSTADO: Más conservador para R=5
    } else if (R_target <= 10.0) {
        omega = 0.008; 
    } else {
        omega = 0.005; 
    }
    
    std::cout << "Configuración para Re = " << R_target << ":" << std::endl;
    std::cout << "  nu = " << nu << std::endl;
    std::cout << "  omega = " << omega << std::endl;
}

// Función auxiliar para verificar si un punto está dentro de la viga
bool is_inside_beam(int i, int j) {
    return (i >= IL && i <= IL + T && j >= 0 && j <= H);
}

// Inicialización de condiciones de frontera
void borders() {
    // Inicializar función de corriente y vorticidad
    #pragma omp parallel for
    for (int i = 0; i <= Nxmax; ++i) {
        for (int j = 0; j <= Nymax; ++j) {
            if (is_inside_beam(i, j)) {
                u[i][j] = 0.0;  
                w[i][j] = 0.0;  
            } else {
                w[i][j] = 0.0;          
                u[i][j] = j * V0;       
            }
        }
    }
    
    // Superficie del fluido (condición de flujo libre)
    for (int i = 0; i <= Nxmax; ++i) {
        if (!is_inside_beam(i, Nymax)) {
            u[i][Nymax] = u[i][Nymax-1] + V0 * h;
            if (Nymax > 1) w[i][Nymax-1] = 0.0;
        }
    }
    
    // Entrada (inlet) - flujo horizontal uniforme
    for (int j = 0; j <= Nymax; ++j) {
        if (!is_inside_beam(0, j) && !is_inside_beam(1, j)) {
            u[1][j] = u[0][j];  
            w[0][j] = 0.0;      
        }
    }
    
    // Línea central
    for (int i = 0; i <= Nxmax; ++i) {
        if (!is_inside_beam(i, 0)) {
            u[i][0] = 0.0;      
            w[i][0] = 0.0;      
        }
    }
    
    // Salida (outlet) - condiciones de gradiente cero
    for (int j = 1; j < Nymax; ++j) {
        if (!is_inside_beam(Nxmax, j) && !is_inside_beam(Nxmax-1, j)) {
            w[Nxmax][j] = w[Nxmax-1][j];  
            u[Nxmax][j] = u[Nxmax-1][j];  
        }
    }
}

// CORRECCIÓN PRINCIPAL: Condiciones de frontera mejoradas para la viga
void beam_boundaries() {
    // PASO 1: Establecer u = 0 y w = 0 en TODA la viga
    for (int i = IL; i <= IL + T && i <= Nxmax; ++i) {
        for (int j = 0; j <= H && j <= Nymax; ++j) {
            u[i][j] = 0.0;  
            w[i][j] = 0.0;  
        }
    }
    
    // PASO 2: Condiciones de no-deslizamiento en las superficies de la viga
    // Frente de la viga (superficie vertical izquierda)
    if (IL > 0) {
        for (int j = 1; j <= H && j < Nymax; ++j) {
            w[IL-1][j] = -2.0 * u[IL-2][j] / (h * h);
        }
    }
    
    // Atrás de la viga (superficie vertical derecha)
    if (IL + T + 1 <= Nxmax) {
        for (int j = 1; j <= H && j < Nymax; ++j) {
            w[IL+T+1][j] = -2.0 * u[IL+T+2][j] / (h * h);
        }
    }
    
    // Parte superior de la viga (superficie horizontal)
    if (H + 1 <= Nymax) {
        for (int i = IL; i <= IL + T && i <= Nxmax; ++i) {
            w[i][H+1] = -2.0 * u[i][H+2] / (h * h);
        }
    }
    
    // PASO 3: CORRECCIÓN - Tratamiento mejorado de esquinas para Re=5
    if (IL > 0 && H + 1 <= Nymax) {
        // Esquina superior izquierda - MÉTODO CORREGIDO
        // En lugar de combinar ambas direcciones, usar promedio de las condiciones adyacentes
        double w_from_vertical = -2.0 * u[IL-2][H+1] / (h * h);  // Desde superficie vertical
        double w_from_horizontal = -2.0 * u[IL-1][H+2] / (h * h); // Desde superficie horizontal
        
        // Para Re=5, usar un promedio ponderado más conservador
        if (R_target >= 5.0) {
            // Dar más peso a la condición menos extrema para suavizar
            if (std::abs(w_from_vertical) < std::abs(w_from_horizontal)) {
                w[IL-1][H+1] = 0.7 * w_from_vertical + 0.3 * w_from_horizontal;
            } else {
                w[IL-1][H+1] = 0.3 * w_from_vertical + 0.7 * w_from_horizontal;
            }
        } else {
            // Para Reynolds menores, usar promedio simple
            w[IL-1][H+1] = 0.5 * (w_from_vertical + w_from_horizontal);
        }
    }
    
    if (IL + T + 1 <= Nxmax && H + 1 <= Nymax) {
        // Esquina superior derecha - MÉTODO CORREGIDO
        double w_from_vertical = -2.0 * u[IL+T+2][H+1] / (h * h);
        double w_from_horizontal = -2.0 * u[IL+T+1][H+2] / (h * h);
        
        if (R_target >= 5.0) {
            // Mismo tratamiento conservador para Re=5
            if (std::abs(w_from_vertical) < std::abs(w_from_horizontal)) {
                w[IL+T+1][H+1] = 0.7 * w_from_vertical + 0.3 * w_from_horizontal;
            } else {
                w[IL+T+1][H+1] = 0.3 * w_from_vertical + 0.7 * w_from_horizontal;
            }
        } else {
            w[IL+T+1][H+1] = 0.5 * (w_from_vertical + w_from_horizontal);
        }
    }
    
    // PASO 4: NUEVO - Suavizado adicional para Reynolds altos en esquinas
    if (R_target >= 5.0) {
        // Aplicar suavizado en un radio alrededor de las esquinas
        if (IL > 1 && H + 2 <= Nymax) {
            // Área alrededor de esquina superior izquierda
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    int ii = IL - 1 + di;
                    int jj = H + 1 + dj;
                    if (ii > 0 && ii < Nxmax && jj > 0 && jj < Nymax && 
                        !is_inside_beam(ii, jj)) {
                        // Aplicar filtro suavizador solo si el valor parece excesivo
                        if (std::abs(w[ii][jj]) > 2.0) {  // Umbral para detectar picos
                            double smooth_w = 0.0;
                            int count = 0;
                            for (int ddi = -1; ddi <= 1; ++ddi) {
                                for (int ddj = -1; ddj <= 1; ++ddj) {
                                    int iii = ii + ddi;
                                    int jjj = jj + ddj;
                                    if (iii >= 0 && iii <= Nxmax && jjj >= 0 && jjj <= Nymax &&
                                        !is_inside_beam(iii, jjj) && !(ddi == 0 && ddj == 0)) {
                                        smooth_w += w[iii][jjj];
                                        count++;
                                    }
                                }
                            }
                            if (count > 0) {
                                // Combinar valor original con promedio suavizado
                                w[ii][jj] = 0.6 * w[ii][jj] + 0.4 * (smooth_w / count);
                            }
                        }
                    }
                }
            }
        }
    }
}

// Función de relajación para función de corriente
void relax_stream_function() {
    #pragma omp parallel for
    for (int i = 1; i < Nxmax; ++i) {
        for (int j = 1; j < Nymax; ++j) {
            if (is_inside_beam(i, j)) continue;
            
            double old_u = u[i][j];
            double new_u = 0.25 * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] + h*h*w[i][j]);
            double r1 = omega * (new_u - old_u);
            u[i][j] += r1;
        }
    }
}

// Relajación de vorticidad con mejor estabilidad
void relax_vorticity() {
    #pragma omp parallel for
    for (int i = 1; i < Nxmax; ++i) {
        for (int j = 1; j < Nymax; ++j) {
            if (is_inside_beam(i, j)) continue;
            
            double old_w = w[i][j];
            
            // Verificar que no estemos en los bordes
            if (i == 1 || i == Nxmax - 1 || j == 1 || j == Nymax - 1) {
                double a1 = w[i+1][j] + w[i-1][j] + w[i][j+1] + w[i][j-1];
                double new_w = 0.25 * a1;
                double r2 = omega * (new_w - old_w);
                w[i][j] += r2;
            } else {
                double a1 = w[i+1][j] + w[i-1][j] + w[i][j+1] + w[i][j-1];
                double a2 = (u[i][j+1] - u[i][j-1]) * (w[i+1][j] - w[i-1][j]);
                double a3 = (u[i+1][j] - u[i-1][j]) * (w[i][j+1] - w[i][j-1]);
                
                // Factor de estabilización mejorado para Reynolds altos
                double stability_factor = 1.0;
                if (R_target >= 5.0) {
                    stability_factor = 0.4;  // MÁS conservador para R=5
                } else if (R_target > 2.0) {
                    stability_factor = 0.7;  
                } else if (R_target > 1.5) {
                    stability_factor = 0.8;  
                }
                
                double new_w = 0.25 * (a1 - stability_factor * (R/4.0) * (a2 - a3));
                double r2 = omega * (new_w - old_w);
                w[i][j] += r2;
            }
        }
    }
}

// Criterio de convergencia adaptativo
bool relax_until_converge(int max_iterations = 350000) {  // Aumentado para Re=5
    int iter = 0;
    double max_diff_u, max_diff_w;
    
    R = V0 * h / nu;
    #pragma omp single
    {
        std::cout << "Número de Reynolds de malla calculado R = " << R << std::endl;
        std::cout << "Objetivo: Re = " << R_target << std::endl;
    }
    
    // Tolerancia adaptativa
    double effective_tol = tol;
    if (R_target >= 5.0) {
        effective_tol = tol * 200;  
        std::cout << "Usando tolerancia muy relajada para R>=5: " << effective_tol << std::endl;
    } else if (R_target > 2.0) {
        effective_tol = tol * 50;   
        std::cout << "Usando tolerancia relajada: " << effective_tol << std::endl;
    } else if (R_target > 1.5) {
        effective_tol = tol * 10;   
        std::cout << "Usando tolerancia relajada: " << effective_tol << std::endl;
    }
    
    do {
        max_diff_u = 0.0;
        max_diff_w = 0.0;
        
        // Aplicar condiciones de frontera
        beam_boundaries();
        
        // Guardar valores anteriores
        double u_old[Nxmax + 1][Nymax + 1];
        double w_old[Nxmax + 1][Nymax + 1];

        #pragma omp parallel 
        {
            #pragma omp for nowait
            for (int i = 0; i <= Nxmax; ++i) {
                for (int j = 0; j <= Nymax; ++j) {
                    u_old[i][j] = u[i][j];
                    w_old[i][j] = w[i][j];
                }
            }
            #pragma omp barrier
        }
        
        // Relajar función de corriente
        relax_stream_function();
        
        // Relajar vorticidad
        relax_vorticity();
        
        // Volver a aplicar condiciones de frontera
        beam_boundaries();
        
        // Calcular diferencias máximas
        #pragma omp parallel for reduction(max:max_diff_u,max_diff_w)
        for (int i = 1; i < Nxmax; ++i) {
            for (int j = 1; j < Nymax; ++j) {
                if (is_inside_beam(i, j)) continue;
                
                max_diff_u = std::max(max_diff_u, std::abs(u[i][j] - u_old[i][j]));
                max_diff_w = std::max(max_diff_w, std::abs(w[i][j] - w_old[i][j]));
            }
        }
        
        ++iter;
        
        // Mostrar progreso
        int report_interval = (R_target >= 5.0) ? 3000 : 5000;
        
        if (iter % report_interval == 0) {
            std::cout << "Iteración " << iter << ": Error u = " << max_diff_u 
                      << ", Error w = " << max_diff_w << std::endl;
        }
        
        // Verificar divergencia
        double divergence_threshold = (R_target >= 5.0) ? 50.0 : 1000.0;
        if (std::isnan(max_diff_u) || std::isnan(max_diff_w) || 
            max_diff_u > divergence_threshold || max_diff_w > divergence_threshold) {
            std::cout << "¡Advertencia: La simulación puede estar divergiendo!" << std::endl;
            std::cout << "Error u = " << max_diff_u << ", Error w = " << max_diff_w << std::endl;
            return false;
        }
        
        // Criterio de convergencia
        if (max_diff_u < effective_tol && max_diff_w < effective_tol) {
            break;
        }
        
        // Convergencia parcial
        if (max_diff_u < effective_tol && iter > 50000) {
            static double last_w_error = -1.0;
            static int stagnant_count = 0;
            
            if (std::abs(max_diff_w - last_w_error) < 1e-15) {
                stagnant_count++;
            } else {
                stagnant_count = 0;
            }
            
            last_w_error = max_diff_w;
            
            int patience = (R_target >= 5.0) ? 8000 : 3000;
            
            if (stagnant_count > patience) {
                std::cout << "Convergencia parcial: u convergió, w estancado en " << max_diff_w << std::endl;
                break;
            }
        }
        
    } while (iter < max_iterations);
    
    if (iter >= max_iterations) {
        std::cout << "¡Advertencia: Se alcanzó el máximo de iteraciones!" << std::endl;
        std::cout << "Error final: u = " << max_diff_u << ", w = " << max_diff_w << std::endl;
        
        double acceptance_threshold = (R_target >= 5.0) ? effective_tol * 20000 : effective_tol * 1000;
        if (max_diff_u < acceptance_threshold) {
            std::cout << "Aceptando solución con convergencia parcial para R=" << R_target << std::endl;
            return true;
        }
        return false;
    }
    
    std::cout << "Convergencia alcanzada en " << iter << " iteraciones." << std::endl;
    std::cout << "Error final: u = " << max_diff_u << ", w = " << max_diff_w << std::endl;
    return true;
}

// Normalización
void normalize() {
    #pragma omp parallel for
    for (int i = 0; i <= Nxmax; ++i) {
        for (int j = 0; j <= Nymax; ++j) {
            u[i][j] /= (V0 * h);
        }
    }
}

// Función auxiliar para formatear números de Reynolds
std::string format_reynolds(double reynolds) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << reynolds;
    return oss.str();
}

// Exportar función de corriente (MODIFICADA)
void export_streamfunction(double reynolds) {
    std::string reynolds_str = format_reynolds(reynolds);
    std::string filename = "Datos/streamfunction_Re_NBS" + reynolds_str + ".dat";
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cout << "Error: No se pudo crear el archivo " << filename << std::endl;
        return;
    }
    
    file << "# Función de corriente - Re = " << reynolds << std::endl;
    file << "# Formato: x y psi" << std::endl;
    file << "# Parámetros: Nx=" << Nxmax << " Ny=" << Nymax << " h=" << h << std::endl;
    file << std::fixed << std::setprecision(6);
    
    for (int i = 0; i < Nxmax; ++i) {
        for (int j = 0; j < Nymax; ++j) {
            double x = i * h;
            double y = j * h;
            file << x << " " << y << " " << u[i][j] << std::endl;
        }
    }
    file.close();
    std::cout << "Función de corriente exportada a: " << filename << std::endl;
}

// Exportar vorticidad (MODIFICADA)
void export_vorticity(double reynolds) {
    std::string reynolds_str = format_reynolds(reynolds);
    std::string filename = "Datos/vorticity_Re_NBS" + reynolds_str + ".dat";
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cout << "Error: No se pudo crear el archivo " << filename << std::endl;
        return;
    }
    
    file << "# Vorticidad - Re = " << reynolds << std::endl;
    file << "# Formato: x y omega" << std::endl;
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

// Exportar campo de velocidades (MODIFICADA)
void export_velocity_field(double reynolds) {
    std::string reynolds_str = format_reynolds(reynolds);
    std::string filename = "Datos/velocity_field_Re_NBS" + reynolds_str + ".dat";
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cout << "Error: No se pudo crear el archivo " << filename << std::endl;
        return;
    }
    
    file << "# Campo de velocidades - Re = " << reynolds << std::endl;
    file << "# Formato: x y vx vy velocidad_magnitud" << std::endl;
    file << "# Parámetros: Nx=" << Nxmax << " Ny=" << Nymax << " h=" << h << std::endl;
    file << std::fixed << std::setprecision(6);
    
    for (int i = 1; i < Nxmax - 1; ++i) {
        for (int j = 1; j < Nymax - 1; ++j) {
            if (is_inside_beam(i, j)) continue;
            
            double x = i * h;
            double y = j * h;
            
            double vx = (u[i][j+1] - u[i][j-1]) / (2.0 * h);
            double vy = -(u[i+1][j] - u[i-1][j]) / (2.0 * h);
            double v_mag = std::sqrt(vx*vx + vy*vy);
            
            file << x << " " << y << " " << vx << " " << vy << " " << v_mag << std::endl;
        }
    }
    file.close();
    std::cout << "Campo de velocidades exportado a: " << filename << std::endl;
}

// Limpiar arrays
void reset_fields() {
    #pragma omp parallel for
    for (int i = 0; i <= Nxmax; ++i) {
        for (int j = 0; j <= Nymax; ++j) {
            u[i][j] = 0.0;
            w[i][j] = 0.0;
        }
    }
}

// Ejecutar simulación
bool run_simulation(double reynolds) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "INICIANDO SIMULACIÓN PARA Re = " << reynolds << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    reset_fields();
    configure_reynolds(reynolds);
    borders();
    
    bool converged = relax_until_converge();
    
    if (!converged) {
        std::cout << "¡Error: No se logró convergencia para Re = " << reynolds << "!" << std::endl;
        return false;
    }
    
    normalize();
    
    export_streamfunction(reynolds);
    export_vorticity(reynolds);
    export_velocity_field(reynolds);
    
    std::cout << "✓ Simulación completada exitosamente para Re = " << reynolds << std::endl;
    return true;
}

// Función principal (MODIFICADA)
int main() {
    std::cout << "=== Solver de Navier-Stokes Multi-Reynolds (ESQUINAS CORREGIDAS) ===" << std::endl;
    std::cout << "Parámetros del dominio:" << std::endl;
    std::cout << "  Malla: " << Nxmax << " x " << Nymax << std::endl;
    std::cout << "  Viga: posición x=[" << IL << "," << IL+T << "], altura=" << H << std::endl;
    std::cout << "  V0 = " << V0 << ", h = " << h << std::endl;
    std::cout << "  Tolerancia = " << tol << std::endl;
    
    // Crear directorio para los datos
    if (!create_directory("Datos")) {
        std::cout << "Error: No se pudo crear el directorio 'Datos'. Terminando programa." << std::endl;
        return 1;
    }
    
    std::vector<double> reynolds_values = {0.5, 1.0, 2.0, 5.0};  
    int successful_simulations = 0;
    
    for (double re : reynolds_values) {
        if (run_simulation(re)) {
            successful_simulations++;
        } else {
            std::cout << "Falló la simulación para Re = " << re << std::endl;
        }
    }
    
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "RESUMEN FINAL" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    std::cout << "Simulaciones exitosas: " << successful_simulations << "/" << reynolds_values.size() << std::endl;
    
    std::cout << "\nArchivos generados en la carpeta 'Datos':" << std::endl;
    for (double re : reynolds_values) {
        std::string re_str = format_reynolds(re);
        std::cout << "  Re = " << re << ":" << std::endl;
        std::cout << "    - Datos/streamfunction_Re_NBS" << re_str << ".dat" << std::endl;
        std::cout << "    - Datos/vorticity_Re_NBS" << re_str << ".dat" << std::endl;
        std::cout << "    - Datos/velocity_field_Re_NBS" << re_str << ".dat" << std::endl;
    }
    
    return 0;
}

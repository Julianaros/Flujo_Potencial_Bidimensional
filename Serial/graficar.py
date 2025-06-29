#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import os
from scipy.interpolate import griddata

# Configuración de matplotlib para mejores gráficas
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'axes.linewidth': 1.2,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10
})

class CFDVisualizer:
    def __init__(self, data_folder="Datos", output_folder="Graficas"):
        """Inicializar el visualizador con parámetros de la simulación"""
        # Parámetros del dominio (deben coincidir con el código C++)
        self.Nxmax = 160
        self.Nymax = 30
        self.IL = 10      # Inicio de la viga
        self.H = 8        # Altura de la viga
        self.T = 8        # Longitud de la viga
        self.h = 1.0      # Espaciado de la malla
        
        # Configuración de carpetas
        self.data_folder = data_folder
        self.output_folder = output_folder
        
        # Crear carpetas si no existen
        self.setup_folders()
        
    def setup_folders(self):
        """Crear carpetas de datos y gráficas si no existen"""
        if not os.path.exists(self.data_folder):
            os.makedirs(self.data_folder)
            print(f"✓ Creada carpeta de datos: {self.data_folder}")
        
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
            print(f"✓ Creada carpeta de gráficas: {self.output_folder}")
    
    def get_data_path(self, filename):
        """Obtener ruta completa del archivo de datos"""
        return os.path.join(self.data_folder, filename)
    
    def get_output_path(self, filename):
        """Obtener ruta completa del archivo de salida"""
        return os.path.join(self.output_folder, filename)
        
    def format_reynolds(self, reynolds):
        """Formatear número de Reynolds para nombres de archivo"""
        return f"{reynolds:.1f}"
        
    def load_data(self, filename):
        """Cargar datos desde archivo .dat en la carpeta de datos"""
        filepath = self.get_data_path(filename)
        try:
            # Leer el archivo saltando las líneas de comentario
            data = pd.read_csv(filepath, sep='\s+', comment='#', header=None)
            if data.shape[1] >= 3:
                x = data.iloc[:, 0].values
                y = data.iloc[:, 1].values
                z = data.iloc[:, 2].values
                return x, y, z
            else:
                print(f"Error: {filepath} no tiene el formato esperado")
                return None, None, None
        except Exception as e:
            print(f"Error al cargar {filepath}: {e}")
            return None, None, None
    
    def load_velocity_data(self, filename):
        """Cargar datos de velocidad (5 columnas: x, y, vx, vy, magnitud)"""
        filepath = self.get_data_path(filename)
        try:
            data = pd.read_csv(filepath, sep='\s+', comment='#', header=None)
            if data.shape[1] >= 5:
                x = data.iloc[:, 0].values
                y = data.iloc[:, 1].values
                vx = data.iloc[:, 2].values
                vy = data.iloc[:, 3].values
                v_mag = data.iloc[:, 4].values
                return x, y, vx, vy, v_mag
            else:
                print(f"Error: {filepath} no tiene el formato esperado para velocidad")
                return None, None, None, None, None
        except Exception as e:
            print(f"Error al cargar {filepath}: {e}")
            return None, None, None, None, None
    
    def create_mesh_grids(self, x, y, z):
        """Crear grillas estructuradas para contornos"""
        # Crear grillas regulares
        xi = np.linspace(x.min(), x.max(), self.Nxmax)
        yi = np.linspace(y.min(), y.max(), self.Nymax)
        Xi, Yi = np.meshgrid(xi, yi)
        
        # Interpolar datos al grid
        Zi = griddata((x, y), z, (Xi, Yi), method='linear', fill_value=0)
        
        return Xi, Yi, Zi
    
    def add_beam_geometry(self, ax):
        """Agregar geometría de la viga al gráfico"""
        beam_rect = patches.Rectangle(
            (self.IL * self.h, 0), 
            self.T * self.h, 
            self.H * self.h,
            linewidth=2, 
            edgecolor='black', 
            facecolor='gray',
            alpha=0.8,
            zorder=10
        )
        ax.add_patch(beam_rect)
        
        # Etiqueta de la viga
        ax.text(
            (self.IL + self.T/2) * self.h, 
            self.H/2 * self.h, 
            'VIGA', 
            ha='center', 
            va='center',
            fontweight='bold',
            fontsize=10,
            color='white',
            zorder=11
        )
    
    def plot_streamlines(self, reynolds, save_fig=True):
        """Graficar función de corriente (líneas de flujo)"""
        re_str = self.format_reynolds(reynolds)
        filename = f"streamfunction_Re{re_str}.dat"
        x, y, psi = self.load_data(filename)
        
        if x is None:
            print(f"No se pudo cargar {filename} desde {self.data_folder}")
            return
        
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Crear grilla
        X, Y, Z = self.create_mesh_grids(x, y, psi)
        
        # Para Re = 5, usar un número de niveles adaptado
        if reynolds >= 5.0:
            num_levels = 25  # Más niveles para capturar complejidad
            contour_levels = 15
        else:
            num_levels = 20
            contour_levels = 12
        
        # Contornos de función de corriente
        levels = np.linspace(np.nanmin(Z), np.nanmax(Z), contour_levels)
        contours = ax.contour(X, Y, Z, levels=levels, colors='blue', linewidths=1.0)
        ax.clabel(contours, inline=True, fontsize=8, fmt='%.2f')
        
        # Contornos rellenos para mejor visualización
        contourf = ax.contourf(X, Y, Z, levels=num_levels, cmap='viridis', alpha=0.6)
        cbar = plt.colorbar(contourf, ax=ax, shrink=0.8)
        cbar.set_label('Función de Corriente ψ', rotation=270, labelpad=20)
        
        # Agregar viga
        self.add_beam_geometry(ax)
        
        # Configuración de ejes
        ax.set_xlim(0, self.Nxmax)
        ax.set_ylim(0, self.Nymax)
        ax.set_xlabel('Posición X')
        ax.set_ylabel('Posición Y')
        ax.set_title(f'Líneas de Flujo (Re = {reynolds})', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Anotar parámetros con nota especial para Re = 5
        textstr = f'Re = {reynolds}\nMalla: {self.Nxmax}×{self.Nymax}'
        if reynolds >= 5.0:
            textstr += '\n(Flujo complejo)'
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        if save_fig:
            output_file = self.get_output_path(f'streamlines_Re{re_str}.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✓ Gráfica de líneas de flujo guardada: {output_file}")
        plt.show()
    
    def plot_vorticity(self, reynolds, save_fig=True):
        """Graficar campo de vorticidad"""
        re_str = self.format_reynolds(reynolds)
        filename = f"vorticity_Re{re_str}.dat"
        x, y, omega = self.load_data(filename)
        
        if x is None:
            print(f"No se pudo cargar {filename} desde {self.data_folder}")
            return
        
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Crear grilla
        X, Y, Z = self.create_mesh_grids(x, y, omega)
        
        # Mapa de colores personalizado para vorticidad (divergente)
        colors = ['#000080', '#4169E1', '#87CEEB', '#FFFFFF', '#FFA07A', '#FF4500', '#8B0000']
        vorticity_cmap = LinearSegmentedColormap.from_list('vorticity', colors, N=256)
        
        # Para Re = 5, ajustar número de niveles y rango
        if reynolds >= 5.0:
            num_levels = 60  # Más niveles para capturar detalles
            contour_lines = 15
        else:
            num_levels = 50
            contour_lines = 10
        
        # Contornos de vorticidad
        v_max = np.nanmax(np.abs(Z))
        levels = np.linspace(-v_max, v_max, num_levels)
        contourf = ax.contourf(X, Y, Z, levels=levels, cmap=vorticity_cmap, extend='both')
        
        # Barra de colores
        cbar = plt.colorbar(contourf, ax=ax, shrink=0.8)
        cbar.set_label('Vorticidad ω [1/s]', rotation=270, labelpad=20)
        
        # Contornos de líneas para mejor definición
        contours = ax.contour(X, Y, Z, levels=contour_lines, colors='black', linewidths=0.5, alpha=0.7)
        
        # Agregar viga
        self.add_beam_geometry(ax)
        
        # Configuración
        ax.set_xlim(0, self.Nxmax)
        ax.set_ylim(0, self.Nymax)
        ax.set_xlabel('Posición X')
        ax.set_ylabel('Posición Y')
        ax.set_title(f'Campo de Vorticidad (Re = {reynolds})', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Anotar parámetros con información adicional para Re = 5
        textstr = f'Re = {reynolds}\nMalla: {self.Nxmax}×{self.Nymax}'
        if reynolds >= 5.0:
            textstr += f'\nVorticidad máx: {v_max:.3f}'
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        if save_fig:
            output_file = self.get_output_path(f'vorticity_Re{re_str}.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✓ Gráfica de vorticidad guardada: {output_file}")
        plt.show()
    
    def plot_velocity_field(self, reynolds, save_fig=True, skip=3):
        """Graficar campo de velocidades - Solo magnitud y líneas de corriente"""
        re_str = self.format_reynolds(reynolds)
        vel_filename = f"velocity_field_Re{re_str}.dat"
        x_vel, y_vel, vx, vy, v_mag = self.load_velocity_data(vel_filename)
        
        if x_vel is None:
            print(f"No se pudo cargar {vel_filename} desde {self.data_folder}")
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
        
        # === GRÁFICA 1: Solo magnitud de velocidad ===
        
        # Crear grilla para magnitud
        X_mag, Y_mag, Z_mag = self.create_mesh_grids(x_vel, y_vel, v_mag)
        
        # Función para crear máscara de la viga
        def create_beam_mask(X, Y):
            """Crear máscara booleana para excluir la viga"""
            beam_mask = np.zeros_like(X, dtype=bool)
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    x_pos = X[i, j] / self.h  # Convertir a coordenadas de malla
                    y_pos = Y[i, j] / self.h
                    # Verificar si está dentro de la viga
                    if (self.IL <= x_pos <= self.IL + self.T and 0 <= y_pos <= self.H):
                        beam_mask[i, j] = True
            return beam_mask
        
        # Aplicar máscara a los datos
        beam_mask = create_beam_mask(X_mag, Y_mag)
        Z_mag_masked = np.where(beam_mask, np.nan, Z_mag)
        
        # Para Re = 5, usar más niveles en los contornos
        if reynolds >= 5.0:
            contour_levels = 60
        else:
            contour_levels = 50
        
        # Contornos de magnitud con colormap plasma
        contourf1 = ax1.contourf(X_mag, Y_mag, Z_mag_masked, levels=contour_levels, cmap='plasma')
        cbar1 = plt.colorbar(contourf1, ax=ax1, shrink=0.8)
        cbar1.set_label('|V| [m/s]', rotation=270, labelpad=20)
        
        # Agregar viga
        self.add_beam_geometry(ax1)
        
        ax1.set_xlim(0, self.Nxmax)
        ax1.set_ylim(0, self.Nymax)
        ax1.set_xlabel('Posición X')
        ax1.set_ylabel('Posición Y')
        ax1.set_title(f'Magnitud de Velocidad (Re = {reynolds})', fontweight='bold')
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        
        # === GRÁFICA 2: Líneas de corriente ===
        
        # Crear grillas regulares para streamplot
        x_unique = np.linspace(x_vel.min(), x_vel.max(), self.Nxmax)
        y_unique = np.linspace(y_vel.min(), y_vel.max(), self.Nymax)
        X_stream, Y_stream = np.meshgrid(x_unique, y_unique)
        
        # Interpolar componentes de velocidad
        VX = griddata((x_vel, y_vel), vx, (X_stream, Y_stream), method='linear', fill_value=0)
        VY = griddata((x_vel, y_vel), vy, (X_stream, Y_stream), method='linear', fill_value=0)
        
        # Aplicar máscara a las componentes de velocidad
        beam_mask_stream = create_beam_mask(X_stream, Y_stream)
        VX_for_stream = np.where(beam_mask_stream, 0, VX)
        VY_for_stream = np.where(beam_mask_stream, 0, VY)
        
        # Crear streamplot
        speed = np.sqrt(VX_for_stream**2 + VY_for_stream**2)
        speed_for_color = np.where(beam_mask_stream, np.nan, speed)
        
        # Ajustar densidad de líneas según Reynolds - MEJORADO para Re = 5
        if reynolds >= 5.0:
            density = 2.0  # Densidad más alta para capturar complejidad
            linewidth = 1.2
            arrowsize = 1.0
        elif reynolds >= 2.0:
            density = 1.8
            linewidth = 1.3
            arrowsize = 1.1
        else:
            density = 1.5
            linewidth = 1.5
            arrowsize = 1.2
        
        strm = ax2.streamplot(X_stream, Y_stream, VX_for_stream, VY_for_stream, 
                            color=speed_for_color, cmap='viridis', 
                            density=density, linewidth=linewidth, arrowsize=arrowsize)
        
        cbar2 = plt.colorbar(strm.lines, ax=ax2, shrink=0.8)
        cbar2.set_label('Velocidad [m/s]', rotation=270, labelpad=20)
        
        # Agregar viga
        self.add_beam_geometry(ax2)
        
        ax2.set_xlim(0, self.Nxmax)
        ax2.set_ylim(0, self.Nymax)
        ax2.set_xlabel('Posición X')
        ax2.set_ylabel('Posición Y')
        ax2.set_title(f'Líneas de Corriente (Re = {reynolds})', fontweight='bold')
        ax2.set_aspect('equal')
        ax2.grid(True, alpha=0.3)
        
        # Anotar parámetros en ambas gráficas
        v_max_val = np.nanmax(Z_mag_masked)
        textstr = f'Re = {reynolds}\nMalla: {self.Nxmax}×{self.Nymax}'
        if reynolds >= 5.0:
            textstr += f'\nV_max: {v_max_val:.3f} m/s'
        
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        ax2.text(0.02, 0.98, textstr, transform=ax2.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        if save_fig:
            output_file = self.get_output_path(f'velocity_field_Re{re_str}.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✓ Gráfica de campo de velocidades guardada: {output_file}")
        plt.show()
    
    def plot_reynolds_comparison(self, reynolds_list, save_fig=True):
        """Crear gráfica comparativa entre diferentes Reynolds"""
        fig, axes = plt.subplots(len(reynolds_list), 3, figsize=(20, 6*len(reynolds_list)))
        
        if len(reynolds_list) == 1:
            axes = axes.reshape(1, -1)
        
        for idx, re in enumerate(reynolds_list):
            re_str = self.format_reynolds(re)
            
            # Cargar datos
            stream_file = f"streamfunction_Re{re_str}.dat"
            vort_file = f"vorticity_Re{re_str}.dat"
            vel_file = f"velocity_field_Re{re_str}.dat"
            
            # Función de corriente
            x_s, y_s, psi = self.load_data(stream_file)
            if x_s is not None:
                X_s, Y_s, Z_s = self.create_mesh_grids(x_s, y_s, psi)
                
                # Ajustar niveles según Reynolds
                if re >= 5.0:
                    num_contours = 20
                else:
                    num_contours = 15
                
                levels = np.linspace(np.nanmin(Z_s), np.nanmax(Z_s), num_contours)
                axes[idx, 0].contour(X_s, Y_s, Z_s, levels=levels, colors='blue', linewidths=0.8)
                contourf_s = axes[idx, 0].contourf(X_s, Y_s, Z_s, levels=30, cmap='viridis', alpha=0.6)
                self.add_beam_geometry(axes[idx, 0])
                axes[idx, 0].set_title(f'Líneas de Flujo - Re = {re}')
                axes[idx, 0].set_aspect('equal')
                
            # Vorticidad
            x_v, y_v, omega = self.load_data(vort_file)
            if x_v is not None:
                X_v, Y_v, Z_v = self.create_mesh_grids(x_v, y_v, omega)
                v_max = np.nanmax(np.abs(Z_v))
                
                # Más niveles para Re = 5
                if re >= 5.0:
                    num_levels = 40
                else:
                    num_levels = 30
                
                levels_v = np.linspace(-v_max, v_max, num_levels)
                contourf_v = axes[idx, 1].contourf(X_v, Y_v, Z_v, levels=levels_v, cmap='RdBu_r')
                self.add_beam_geometry(axes[idx, 1])
                axes[idx, 1].set_title(f'Vorticidad - Re = {re}')
                axes[idx, 1].set_aspect('equal')
                
            # Campo de velocidades
            x_vel, y_vel, vx, vy, v_mag = self.load_velocity_data(vel_file)
            if x_vel is not None:
                X_mag, Y_mag, Z_mag = self.create_mesh_grids(x_vel, y_vel, v_mag)
                # Crear máscara para la viga
                beam_mask = np.zeros_like(X_mag, dtype=bool)
                for i in range(X_mag.shape[0]):
                    for j in range(X_mag.shape[1]):
                        x_pos = X_mag[i, j] / self.h
                        y_pos = Y_mag[i, j] / self.h
                        if (self.IL <= x_pos <= self.IL + self.T and 0 <= y_pos <= self.H):
                            beam_mask[i, j] = True
                
                Z_mag_masked = np.where(beam_mask, np.nan, Z_mag)
                
                # Más niveles para Re = 5
                if re >= 5.0:
                    contour_levels = 40
                else:
                    contour_levels = 30
                    
                contourf_mag = axes[idx, 2].contourf(X_mag, Y_mag, Z_mag_masked, 
                                                   levels=contour_levels, cmap='plasma')
                self.add_beam_geometry(axes[idx, 2])
                axes[idx, 2].set_title(f'Magnitud Velocidad - Re = {re}')
                axes[idx, 2].set_aspect('equal')
            
            # Configurar ejes
            for ax in axes[idx, :]:
                ax.set_xlim(0, self.Nxmax)
                ax.set_ylim(0, self.Nymax)
                ax.grid(True, alpha=0.3)
                if idx == len(reynolds_list) - 1:  # Solo en la última fila
                    ax.set_xlabel('Posición X')
                if ax == axes[idx, 0]:  # Solo en la primera columna
                    ax.set_ylabel('Posición Y')
        
        plt.tight_layout()
        if save_fig:
            re_str = "_".join([self.format_reynolds(re) for re in reynolds_list])
            output_file = self.get_output_path(f'reynolds_comparison_{re_str}.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✓ Gráfica comparativa guardada: {output_file}")
        plt.show()
    
    def plot_all_reynolds(self, reynolds_list=[0.5, 1.0, 2.0, 5.0], save_figs=True):
        """Graficar todos los casos de Reynolds incluyendo Re = 5.0"""
        print("=" * 60)
        print("  VISUALIZADOR DE RESULTADOS CFD")
        print("  Simulación de Flujo Alrededor de Viga")
        print("  ORGANIZADO CON CARPETAS")
        print("=" * 60)
        print(f"  Carpeta de datos: {self.data_folder}")
        print(f"  Carpeta de gráficas: {self.output_folder}")
        print("=" * 60)
        
        valid_reynolds = []
        
        for re in reynolds_list:
            print(f"\nProcesando Reynolds = {re}")
            re_str = self.format_reynolds(re)
            
            # Verificar que existan los archivos
            files_to_check = [
                f"streamfunction_Re{re_str}.dat",
                f"vorticity_Re{re_str}.dat", 
                f"velocity_field_Re{re_str}.dat"
            ]
            
            files_exist = [os.path.exists(self.get_data_path(f)) for f in files_to_check]
            
            if not all(files_exist):
                missing_files = [f for f, exists in zip(files_to_check, files_exist) if not exists]
                print(f"Faltan archivos para Re = {re}: {missing_files}")
                print(f"  Buscando en: {self.data_folder}")
                if re == 5.0:
                    print("  NOTA: Re = 5.0 requiere convergencia especial. Verifique que el solver C++ haya completado exitosamente.")
                continue
            
            try:
                print("  → Graficando líneas de flujo...")
                self.plot_streamlines(re, save_figs)
                
                print("  → Graficando campo de velocidades...")
                self.plot_velocity_field(re, save_figs)
                
                print("  → Graficando vorticidad...")
                self.plot_vorticity(re, save_figs)
                
                valid_reynolds.append(re)
                print(f"✓ Completado Re = {re}")
                
                # Información adicional para Re = 5
                if re == 5.0:
                    print("  → Re = 5.0: Flujo con mayor complejidad, revise patrones de recirculación")
                
            except Exception as e:
                print(f"✗ Error procesando Re = {re}: {e}")
                if re == 5.0:
                    print("  → Para Re = 5.0, verifique la convergencia del solver y la calidad de los datos")
        
        # Crear gráfica comparativa si hay múltiples Reynolds válidos
        if len(valid_reynolds) > 1:
            print(f"\nCreando gráfica comparativa para Re = {valid_reynolds}")
            try:
                self.plot_reynolds_comparison(valid_reynolds, save_figs)
                print("✓ Gráfica comparativa completada")
            except Exception as e:
                print(f"✗ Error en gráfica comparativa: {e}")
        
        # Resumen final
        print(f"\n{'='*60}")
        print("RESUMEN DE VISUALIZACIÓN")
        print(f"{'='*60}")
        print(f"Reynolds procesados exitosamente: {len(valid_reynolds)}/{len(reynolds_list)}")
        print(f"Datos leídos desde: {os.path.abspath(self.data_folder)}")
        print(f"Gráficas guardadas en: {os.path.abspath(self.output_folder)}")
        

def main():
    """Función principal"""
    # Crear instancia del visualizador con carpetas organizadas
    visualizer = CFDVisualizer(data_folder="Datos", output_folder="Graficas")
    
    # Lista de valores de Reynolds - AMPLIADA para incluir Re = 5.0
    reynolds_values = [0.5, 1.0, 2.0, 5.0]
    
    print("Iniciando visualización con organización de carpetas")
    print(f"Reynolds a procesar: {reynolds_values}")
    print("Estructura de carpetas:")
    print("  Datos/     <- Archivos .dat de entrada")
    print("  Graficas/  <- Archivos .png de salida")
    
    # Generar todas las visualizaciones
    visualizer.plot_all_reynolds(reynolds_values, save_figs=True)

if __name__ == "__main__":
    main()
    
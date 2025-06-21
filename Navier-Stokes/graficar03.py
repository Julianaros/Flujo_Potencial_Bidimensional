#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import os

# Configuración de matplotlib para mejores gráficas
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'axes.linewidth': 1.2,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

class CFDVisualizer:
    def __init__(self):
        """Inicializar el visualizador con parámetros de la simulación"""
        # Parámetros del dominio (deben coincidir con el código C++)
        self.Nxmax = 70
        self.Nymax = 20
        self.IL = 10      # Inicio de la viga
        self.H = 8        # Altura de la viga
        self.T = 8        # Longitud de la viga
        self.h = 1.0      # Espaciado de la malla
        
        # Datos
        self.streamfunction_data = None
        self.vorticity_data = None
        self.velocity_data = None
        
    def load_data(self):
        """Cargar datos de los archivos generados por la simulación"""
        try:
            # Cargar función de corriente
            if os.path.exists('streamfunction.dat'):
                self.streamfunction_data = pd.read_csv('streamfunction.dat', 
                                                     sep=' ', 
                                                     comment='#',
                                                     names=['x', 'y', 'psi'])
                print(f"✓ Función de corriente cargada: {len(self.streamfunction_data)} puntos")
            
            # Cargar vorticidad
            if os.path.exists('vorticity.dat'):
                self.vorticity_data = pd.read_csv('vorticity.dat', 
                                                sep=' ', 
                                                comment='#',
                                                names=['x', 'y', 'omega'])
                print(f"✓ Vorticidad cargada: {len(self.vorticity_data)} puntos")
            
            # Cargar campo de velocidades
            if os.path.exists('velocity_field.dat'):
                self.velocity_data = pd.read_csv('velocity_field.dat', 
                                               sep=' ', 
                                               comment='#',
                                               names=['x', 'y', 'vx', 'vy', 'v_mag'])
                print(f"✓ Campo de velocidades cargado: {len(self.velocity_data)} puntos")
                
        except Exception as e:
            print(f"Error al cargar datos: {e}")
            print("Asegúrate de que los archivos .dat estén en el directorio actual")
    
    def create_mesh_grids(self, data):
        """Crear grillas para contornos a partir de datos dispersos"""
        x_unique = np.sort(data['x'].unique())
        y_unique = np.sort(data['y'].unique())
        X, Y = np.meshgrid(x_unique, y_unique)
        
        # Crear matriz Z
        Z = np.full(X.shape, np.nan)
        for _, row in data.iterrows():
            i = np.where(y_unique == row['y'])[0][0]
            j = np.where(x_unique == row['x'])[0][0]
            Z[i, j] = row.iloc[2]  # Tercera columna (psi, omega, etc.)
        
        return X, Y, Z
    
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
    
    def plot_streamfunction(self, save_fig=True):
        """Graficar función de corriente (líneas de flujo)"""
        if self.streamfunction_data is None:
            print(" No hay datos de función de corriente")
            return
        
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Crear grilla
        X, Y, Z = self.create_mesh_grids(self.streamfunction_data)
        
        # Contornos de función de corriente
        levels = np.linspace(np.nanmin(Z), np.nanmax(Z), 20)
        contours = ax.contour(X, Y, Z, levels=levels, colors='blue', linewidths=1.0)
        ax.clabel(contours, inline=True, fontsize=8, fmt='%.2f')
        
        # Contornos rellenos para mejor visualización
        contourf = ax.contourf(X, Y, Z, levels=50, cmap='viridis', alpha=0.6)
        cbar = plt.colorbar(contourf, ax=ax, shrink=0.8)
        cbar.set_label('Función de Corriente ψ', rotation=270, labelpad=20)
        
        # Agregar viga
        self.add_beam_geometry(ax)
        
        # Configuración de ejes
        ax.set_xlabel('Posición X')
        ax.set_ylabel('Posición Y')
        ax.set_title('Líneas de Flujo (Función de Corriente)', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        # Anotar parámetros
        textstr = f'Re_malla = {1.0 * self.h / 1.0:.1f}\nMalla: {self.Nxmax}×{self.Nymax}'
        props = dict(boxstyle='round', facecolor='white', alpha=0.8)
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        if save_fig:
            plt.savefig('streamlines.png', dpi=300, bbox_inches='tight')
            print("✓ Gráfica de líneas de flujo guardada: streamlines.png")
        plt.show()
    
    def plot_vorticity(self, save_fig=True):
        """Graficar campo de vorticidad"""
        if self.vorticity_data is None:
            print(" No hay datos de vorticidad")
            return
        
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Crear grilla
        X, Y, Z = self.create_mesh_grids(self.vorticity_data)
        
        # Mapa de colores personalizado para vorticidad
        colors = ['blue', 'cyan', 'white', 'yellow', 'red']
        n_bins = 100
        cmap = LinearSegmentedColormap.from_list('vorticity', colors, N=n_bins)
        
        # Contornos de vorticidad
        v_max = np.nanmax(np.abs(Z))
        levels = np.linspace(-v_max, v_max, 50)
        contourf = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, extend='both')
        
        # Barra de colores
        cbar = plt.colorbar(contourf, ax=ax, shrink=0.8)
        cbar.set_label('Vorticidad ω [1/s]', rotation=270, labelpad=20)
        
        # Contornos de líneas
        contours = ax.contour(X, Y, Z, levels=10, colors='black', linewidths=0.5, alpha=0.7)
        
        # Agregar viga
        self.add_beam_geometry(ax)
        
        # Configuración
        ax.set_xlabel('Posición X')
        ax.set_ylabel('Posición Y')
        ax.set_title('Campo de Vorticidad', fontsize=14, fontweight='bold')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        if save_fig:
            plt.savefig('vorticity.png', dpi=300, bbox_inches='tight')
            print("✓ Gráfica de vorticidad guardada: vorticity.png")
        plt.show()
    
    def plot_velocity_field(self, save_fig=True, skip=3):
        """Graficar campo de velocidades con vectores"""
        if self.velocity_data is None:
            print(" No hay datos de campo de velocidades")
            return
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
        
        # --- Gráfica 1: Magnitud de velocidad con vectores ---
        
        # Crear grilla para magnitud
        data_mag = self.velocity_data[['x', 'y', 'v_mag']].copy()
        X, Y, Z_mag = self.create_mesh_grids(data_mag)
        
        # Contornos de magnitud
        contourf1 = ax1.contourf(X, Y, Z_mag, levels=50, cmap='plasma')
        cbar1 = plt.colorbar(contourf1, ax=ax1, shrink=0.8)
        cbar1.set_label('|V| [m/s]', rotation=270, labelpad=20)
        
        # Vectores de velocidad (subsampling para claridad)
        vel_sub = self.velocity_data.iloc[::skip*skip]  # Reducir densidad de vectores
        scale_factor = 30  # Ajustar según sea necesario
        
        ax1.quiver(vel_sub['x'], vel_sub['y'], 
                  vel_sub['vx'], vel_sub['vy'],
                  scale=scale_factor, scale_units='xy', angles='xy',
                  color='white', alpha=0.7, width=0.003)
        
        # Agregar viga
        self.add_beam_geometry(ax1)
        
        ax1.set_xlabel('Posición X')
        ax1.set_ylabel('Posición Y')
        ax1.set_title('Magnitud de Velocidad + Vectores', fontweight='bold')
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        
        # --- Gráfica 2: Líneas de corriente con magnitud ---
        
        # Streamplot usando componentes de velocidad
        x_stream = self.velocity_data['x'].values
        y_stream = self.velocity_data['y'].values
        vx_stream = self.velocity_data['vx'].values
        vy_stream = self.velocity_data['vy'].values
        
        # Crear grilla regular para streamplot
        x_unique = np.sort(self.velocity_data['x'].unique())
        y_unique = np.sort(self.velocity_data['y'].unique())
        X_stream, Y_stream = np.meshgrid(x_unique, y_unique)
        
        VX = np.full(X_stream.shape, np.nan)
        VY = np.full(X_stream.shape, np.nan)
        
        for _, row in self.velocity_data.iterrows():
            i = np.where(y_unique == row['y'])[0][0]
            j = np.where(x_unique == row['x'])[0][0]
            VX[i, j] = row['vx']
            VY[i, j] = row['vy']
        
        # Rellenar NaN con interpolación
        mask = ~np.isnan(VX)
        if np.any(mask):
            # Líneas de corriente
            speed = np.sqrt(VX**2 + VY**2)
            strm = ax2.streamplot(X_stream, Y_stream, VX, VY, 
                                color=speed, cmap='viridis', 
                                density=2, linewidth=1.5)
            
            cbar2 = plt.colorbar(strm.lines, ax=ax2, shrink=0.8)
            cbar2.set_label('Velocidad [m/s]', rotation=270, labelpad=20)
        
        # Agregar viga
        self.add_beam_geometry(ax2)
        
        ax2.set_xlabel('Posición X')
        ax2.set_ylabel('Posición Y')
        ax2.set_title('Líneas de Corriente', fontweight='bold')
        ax2.set_aspect('equal')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        if save_fig:
            plt.savefig('velocity_field.png', dpi=300, bbox_inches='tight')
            print("✓ Gráfica de campo de velocidades guardada: velocity_field.png")
        plt.show()
    


def main():
    """Función principal"""
    print("=" * 60)
    print("  VISUALIZADOR DE RESULTADOS CFD")
    print("  Simulación de Flujo Alrededor de Viga")
    print("=" * 60)
    
    # Crear visualizador
    viz = CFDVisualizer()
    
    # Cargar datos
    print("\n1. Cargando datos...")
    viz.load_data()
    
    # Crear visualizaciones
    print("\n2. Generando visualizaciones...")
    
    # Gráficas individuales
    viz.plot_streamfunction()
    viz.plot_vorticity()
    viz.plot_velocity_field()
    
    print("\n¡Visualización completada exitosamente!")
    print("\nArchivos generados:")
    print("• streamlines.png - Líneas de flujo")
    print("• vorticity.png - Campo de vorticidad")  
    print("• velocity_field.png - Campo de velocidades")

if __name__ == "__main__":
    main()
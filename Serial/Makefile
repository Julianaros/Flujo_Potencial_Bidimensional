# Compilar el archivo beam.cpp
beam:
	@echo "Compilando beam.cpp..."
	g++ beam.cpp -o beam

# Ejecutar el archivo binario
run: beam
	@echo "Ejecutando el binario ./beam..."
	time ./beam

# Ejecutar el script de Python
plot:
	@echo "Ejecutando el script graficar.py..."
	python graficar.py

# Eliminar carpetas Datos y Graficas si existen
clean_dirs:
	@echo "Eliminando carpetas Datos y Graficas..."
	rm -rf Datos Graficas

# Eliminar binarios
clean_bin:
	@echo "Eliminando binario beam..."
	rm -f beam

# Limpiar todo: carpetas + binarios
clean: clean_dirs clean_bin
	@echo "Limpieza completa."

# Ejecutar todo el proceso: compilar, correr y graficar
all: run plot
	@echo "Proceso completo."

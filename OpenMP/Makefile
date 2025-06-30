# Compilador y flags
CXX = g++
CXXFLAGS = -fopenmp -O2

# Archivos fuente y binarios
SOURCES = beam-collapse.cpp beam-dynamic.cpp beam-parallel-for.cpp beam-static.cpp
BINARIES = beam-collapse beam-dynamic beam-parallel-for beam-static

# Compilación individual
beam-collapse: beam-collapse.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

beam-dynamic: beam-dynamic.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

beam-parallel-for: beam-parallel-for.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

beam-static: beam-static.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Compilar todos
build: $(BINARIES)
	@echo "Todos los binarios han sido compilados."

# Ejecutar todos
run: build
	@echo "Ejecutando programas..."
	@bash -c "time ./beam-collapse"
	@bash -c "time ./beam-dynamic"
	@bash -c "time ./beam-parallel-for"
	@bash -c "time ./beam-static"

# Graficar
plot:
	@echo "Ejecutando scripts de graficación..."
	python graficar-collapse.py
	python graficar-dynamic.py
	python graficar_parallelfor.py
	python graficar-static.py

# Limpiar carpetas
clean_dirs:
	@echo "Eliminando carpetas Datos y Graficas..."
	rm -rf Datos Graficas

# Limpiar binarios
clean_bin:
	@echo "Eliminando binarios..."
	rm -f $(BINARIES)

# Limpiar todo
clean: clean_dirs clean_bin
	@echo "Limpieza completa."

# Ejecutar todo
all: run plot
	@echo "Proceso completo."

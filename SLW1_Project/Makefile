# Compilador e flags
FC = gfortran
FFLAGS = -O2 -Wall -Jobj -Iobj

# Lista de módulos e programa principal
MODULES = src/modules/slw1_io.f95 \
          src/modules/slw1_physics.f95 \
          src/modules/radiation_solver.f95

MAIN = src/main.f95

# Objetos correspondentes
OBJS = $(patsubst src/%.f95, obj/%.o, $(MODULES)) obj/main.o

# Regra padrão
all: bin/slw1.exe

# Compilação do executável
bin/slw1.exe: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# Compilar os módulos
obj/%.o: src/%.f95
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c -o $@ $<

# Compilar o main
obj/main.o: $(MAIN)
	$(FC) $(FFLAGS) -c -o $@ $<

# Limpeza
clean:
	rm -rf obj/*.o obj/*.mod bin/slw1.exe

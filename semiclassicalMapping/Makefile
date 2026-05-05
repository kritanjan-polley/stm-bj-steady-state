FC = gfortran
CC = gcc

TARGET = main

FFLAGS = -O2 -ffast-math -funroll-loops -cpp
CFLAGS = -O2
LDFLAGS = -framework Accelerate

# Fortran source codes in order
FSRC = params.f90 \
        functions.f90 \
        propagator.f90 \
        initial.f90 \
        main.f90

# C sources
CSRC = ../lindbladMasterEq/memory.c

OBJ = $(FSRC:.f90=.o) $(notdir $(CSRC:.c=.o))

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: ../lindbladMasterEq/%.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod $(TARGET)

.PHONY: all clean

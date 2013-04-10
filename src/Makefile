CC=gcc
CFLAGS=-g -Wall

all: harmonic_verlet_anderson

harmonic_verlet_anderson: *.c *.h
	$(CC) $(CFLAGS) -c -DHARMONIC -DANDERSON -c -o util.o util.c
	$(CC) $(CFLAGS) -c -DANDERSON -o thermostat.o anderson_thermostat.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DHARMONIC -o force.o harmonic_force.c
	$(CC) $(CFLAGS) -o harmonic_verlet_anderson thermostat.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv harmonic_verlet_anderson ..

clean: 
	rm -rf *o 
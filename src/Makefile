CC=gcc
CFLAGS=-g -Wall -O3

all: harmonic_verlet_anderson lj_verlet_anderson lj_verlet

min_image: min_image.c min_image.h
	$(CC) $(CFLAGS) -c -o min_image.o min_image.c

anderson: util.c util.h anderson_thermostat.c thermostat.h 
	$(CC) $(CFLAGS) -c -DANDERSON -o thermostat.o anderson_thermostat.c

bussi: util.c util.h bussi_thermostat.c thermostat.h 
	$(CC) $(CFLAGS) -c -DBUSSI -o thermostat.o bussi_thermostat.c


harmonic_verlet_anderson: *.c *.h min_image anderson
	$(CC) $(CFLAGS) -c -DHARMONIC -DANDERSON -o util.o util.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DHARMONIC -o force.o harmonic_force.c
	$(CC) $(CFLAGS) -DANDERSON -o harmonic_verlet_anderson min_image.o thermostat.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv harmonic_verlet_anderson ..

lj_verlet_bussi: *.c *.h min_image bussi
	$(CC) $(CFLAGS) -c -DLJ -DBUSSI -o util.o util.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DLJ -o force.o lj_force.c
	$(CC) $(CFLAGS) -DBUSSI -o lj_verlet_bussi min_image.o thermostat.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv lj_verlet_bussi ..


lj_verlet_anderson: *.c *.h min_image anderson
	$(CC) $(CFLAGS) -c -DLJ -DANDERSON -o util.o util.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DLJ -o force.o lj_force.c
	$(CC) $(CFLAGS) -DANDERSON -o lj_verlet_anderson min_image.o thermostat.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv lj_verlet_anderson ..

lj_verlet: *.c *.h min_image
	$(CC) $(CFLAGS) -c -DLJ -o util.o util.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DLJ -o force.o lj_force.c
	$(CC) $(CFLAGS) -o lj_verlet min_image.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv lj_verlet ..

clean: 
	rm -rf *o

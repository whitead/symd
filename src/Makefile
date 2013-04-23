CC=gcc
CFLAGS=-g -Wall -O3

all: harmonic_vverlet_anderson lj_vverlet_anderson lj_vverlet lj_vverlet_bussi harmonic_vverlet_bussi harmonic_vverlet

min_image: min_image.c min_image.h
	$(CC) $(CFLAGS) -c -o min_image.o min_image.c

min_image_nopbc: min_image.c min_image.h
	$(CC) $(CFLAGS) -DNO_PBC -c -o min_image_nopbc.o min_image.c

vverlet: vverlet_integrate.c integrate.h
	$(CC) $(CFLAGS) -c -o vverlet_integrate.o vverlet_integrate.c

harmonic_force: harmonic_force.c force.h
	$(CC) $(CFLAGS) -c -DHARMONIC -o harmonic_force.o harmonic_force.c

anderson: util.c util.h anderson_thermostat.c thermostat.h 
	$(CC) $(CFLAGS) -c -DANDERSON -o anderson_thermostat.o anderson_thermostat.c

bussi: util.c util.h bussi_thermostat.c thermostat.h 
	$(CC) $(CFLAGS) -c -DBUSSI -o bussi_thermostat.o bussi_thermostat.c

harmonic_vverlet: *.c *.h min_image_nopbc vverlet harmonic_force
	$(CC) $(CFLAGS) -c -DHARMONIC -DNO_PBC -o util_harmonic_nopbc.o util.c
	$(CC) $(CFLAGS) -o harmonic_vverlet min_image_nopbc.o util_harmonic_nopbc.o harmonic_force.o vverlet_integrate.o -lgsl -lblas -lm main_loop.c
	mv harmonic_vverlet ..

harmonic_vverlet_anderson: *.c *.h min_image_nopbc anderson vverlet harmonic_force
	$(CC) $(CFLAGS) -c -DHARMONIC -DNO_PBC -DANDERSON -o util_harmonic_nopbc_anderson.o util.c
	$(CC) $(CFLAGS) -DANDERSON -o harmonic_vverlet_anderson min_image_nopbc.o anderson_thermostat.o util_harmonic_nopbc_anderson.o harmonic_force.o vverlet_integrate.o -lgsl -lblas -lm main_loop.c
	mv harmonic_vverlet_anderson ..


lj_vverlet_bussi: *.c *.h min_image bussi
	$(CC) $(CFLAGS) -c -DLJ -DBUSSI -o util.o util.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DLJ -o force.o lj_force.c
	$(CC) $(CFLAGS) -DBUSSI -o lj_vverlet_bussi min_image.o thermostat.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv lj_vverlet_bussi ..

harmonic_vverlet_bussi: *.c *.h min_image_nopbc bussi vverlet harmonic_force
	$(CC) $(CFLAGS) -c -DHARMONIC -DNO_PBC -DBUSSI -o util_harmonic_nopbc_bussi.o util.c
	$(CC) $(CFLAGS) -DBUSSI -o harmonic_vverlet_bussi min_image_nopbc.o bussi_thermostat.o util_harmonic_nopbc_bussi.o harmonic_force.o vverlet_integrate.o -lgsl -lblas -lm main_loop.c
	mv harmonic_vverlet_bussi ..


lj_vverlet_anderson: *.c *.h min_image anderson
	$(CC) $(CFLAGS) -c -DLJ -DANDERSON -o util.o util.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DLJ -o force.o lj_force.c
	$(CC) $(CFLAGS) -DANDERSON -o lj_vverlet_anderson min_image.o thermostat.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv lj_vverlet_anderson ..

lj_vverlet: *.c *.h min_image
	$(CC) $(CFLAGS) -c -DLJ -o util.o util.c
	$(CC) $(CFLAGS) -c -o integrate.o vverlet_integrate.c
	$(CC) $(CFLAGS) -c -DLJ -o force.o lj_force.c
	$(CC) $(CFLAGS) -o lj_vverlet min_image.o util.o force.o integrate.o -lgsl -lblas -lm main_loop.c
	mv lj_vverlet ..

clean: 
	rm -rf *o

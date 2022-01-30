CC=gcc
CFLAGS=-g -Wall -Wextra -Wno-unknown-pragmas -Wno-unused-parameter
PPFLAGS=
SUB_CFLAGS= $(CFLAGS) $(PPFLAGS) -o $@_$(SUFFIX).o -c
FINISH=mv $@ ..
MAKE_EXE_FLAGS=
MAKE_EXE=$(CC) $(CFLAGS) $(PPFLAGS) $(MAKE_EXE_FLAGS) -o $@ $(UTIL_O) $(addsuffix _$(SUFFIX).o, $^) -lgsl -lgslcblas -lm main_loop.c


.PHONY: clean all all-single all all-debug

all-single:
	make optimize harmonic_vverlet
	make optimize lj_vverlet
	make optimize harmonic_vverlet_anderson
	make optimize lj_vverlet_bussi
	make optimize lj_vverlet_anderson
	make optimize harmonic_vverlet_bussi
	make optimize soft_vverlet

all:
	make openmp optimize harmonic_vverlet
	make openmp optimize lj_vverlet
	make openmp optimize harmonic_vverlet_anderson
	make openmp optimize lj_vverlet_bussi
	make openmp optimize lj_vverlet_anderson
	make openmp optimize harmonic_vverlet_bussi
	make openmp optimize soft_vverlet
	make clean

all-debug:
	make harmonic_vverlet
	make lj_vverlet
	make harmonic_vverlet_anderson
	make lj_vverlet_bussi
	make lj_vverlet_anderson
	make harmonic_vverlet_bussi
	make soft_vverlet


cJSON: cJSON.c cJSON.h
	$(CC) $(SUB_CFLAGS) cJSON.c

nlist: nlist.c nlist.h
	$(CC) $(SUB_CFLAGS) nlist.c

min_image: min_image.c min_image.h
	$(CC) $(SUB_CFLAGS) min_image.c

vverlet: vverlet_integrate.c integrate.h
	$(CC) $(SUB_CFLAGS) vverlet_integrate.c

forces: harmonic_force.c lj_forces.c force.h
	$(CC) $(SUB_CFLAGS) harmonic_force.c

soft_force: soft_force.c force.h
	$(CC) $(SUB_CFLAGS) soft_force.c

anderson_thermostat: util.c util.h anderson_thermostat.c thermostat.h
	$(CC) $(SUB_CFLAGS) -DANDERSON anderson_thermostat.c

bussi_thermostat: util.c util.h bussi_thermostat.c thermostat.h
	$(CC) $(SUB_CFLAGS) -DBUSSI bussi_thermostat.c

util: util.c util.h *.h
	$(CC) $(SUB_CFLAGS) util.c

openmp: CFLAGS += -fopenmp
openmp: $(filter-out openmp, $(MAKECMDGOALS))


optimize: CFLAGS += -march=native -Ofast -fmodulo-sched-allow-regmoves
optimize: MAKE_EXE_FLAGS += -fwhole-program
optimize: $(filter-out optimize, $(MAKECMDGOALS))

soft_vverlet: SUFFIX = s
soft_vverlet: min_image soft_force vverlet util cJSON
	$(MAKE_EXE)
	$(FINISH)

harmonic_vverlet: PPFLAGS += -DHARMONIC
harmonic_vverlet: SUFFIX = hv
harmonic_vverlet: min_image vverlet forces nlist util cJSON
	$(MAKE_EXE)
	$(FINISH)


lj_vverlet: PPFLAGS += -DLJ
lj_vverlet: SUFFIX = ljv
lj_vverlet: min_image forces vverlet nlist util cJSON
	$(MAKE_EXE)
	$(FINISH)

harmonic_vverlet_anderson: PPFLAGS += -DHARMONIC -DANDERSON
harmonic_vverlet_anderson: SUFFIX = hva
harmonic_vverlet_anderson: min_image anderson_thermostat vverlet forces nlist util cJSON
	$(MAKE_EXE)
	$(FINISH)

lj_vverlet_bussi: PPFLAGS += -DLJ -DBUSSI
lj_vverlet_bussi: SUFFIX = ljvb
lj_vverlet_bussi: min_image bussi_thermostat forces vverlet nlist util cJSON
	$(MAKE_EXE)
	$(FINISH)


lj_vverlet_anderson: PPFLAGS += -DLJ -DANDERSON
lj_vverlet_anderson: SUFFIX = ljvb
lj_vverlet_anderson: min_image anderson_thermostat forces vverlet nlist util cJSON
	$(MAKE_EXE)
	$(FINISH)


harmonic_vverlet_bussi: PPFLAGS += -DHARMONIC -DBUSSI
harmonic_vverlet_bussi: SUFFIX = hvb
harmonic_vverlet_bussi: min_image bussi_thermostat forces vverlet util cJSON nlist
	$(MAKE_EXE)
	$(FINISH)


clean:
	rm -f *.o


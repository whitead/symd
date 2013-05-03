CC=gcc
CFLAGS=-g -Wall
SUB_CFLAGS= $(CFLAGS) -o $@_$(SUFFIX).o -c
FINISH=mv $@ ..
MAKE_EXE=$(CC) $(CFLAGS) -o $@ $(UTIL_O) $(addsuffix _$(SUFFIX).o, $^) -lgsl -lblas -lm main_loop.c


.PHONY: clean all all-profile all-optimize

all:
	make harmonic_vverlet
	make lj_vverlet
	make harmonic_vverlet_anderson
	make lj_vverlet_bussi
	make lj_vverlet_anderson
	make harmonic_vverlet_bussi

all-profile:
	make profile harmonic_vverlet
	make profile lj_vverlet
	make profile harmonic_vverlet_anderson
	make profile lj_vverlet_bussi
	make profile lj_vverlet_anderson
	make profile harmonic_vverlet_bussi

all-optimize:
	make optimize harmonic_vverlet
	make optimize lj_vverlet
	make optimize harmonic_vverlet_anderson
	make optimize lj_vverlet_bussi
	make optimize lj_vverlet_anderson
	make optimize harmonic_vverlet_bussi



nlist: nlist.c nlist.h
	$(CC) $(SUB_CFLAGS) nlist.c

min_image: min_image.c min_image.h
	$(CC) $(SUB_CFLAGS) min_image.c

min_image_nopbc: min_image.c min_image.h
	$(CC) $(SUB_CFLAGS) -DNO_PBC min_image.c

vverlet: vverlet_integrate.c integrate.h
	$(CC) $(SUB_CFLAGS) vverlet_integrate.c

harmonic_force: harmonic_force.c force.h
	$(CC) $(SUB_CFLAGS) -DHARMONIC harmonic_force.c


lj_force: lj_force.c force.h
	$(CC) $(SUB_CFLAGS) -DLJ lj_force.c

anderson_thermostat: util.c util.h anderson_thermostat.c thermostat.h 
	$(CC) $(SUB_CFLAGS) -DANDERSON anderson_thermostat.c

bussi_thermostat: util.c util.h bussi_thermostat.c thermostat.h 
	$(CC) $(SUB_CFLAGS) -DBUSSI bussi_thermostat.c

util: util.c util.h *.h
	$(CC) $(SUB_CFLAGS) util.c


profile: CFLAGS += -fprofile-generate -march=native -Ofast
profile: FINISH = cd profile && ../$@ run_params_$(PARAMS).txt
profile: $(filter-out profile, $(MAKECMDGOALS))

optimize: CFLAGS += -fprofile-use -march=native -Ofast
optimize: $(filter-out optimize, $(MAKECMDGOALS))

harmonic_vverlet: SUB_CFLAGS += -DHARMONIC -DNO_PBC
harmonic_vverlet: SUFFIX = hv
harmonic_vverlet: PARAMS=harmonic
harmonic_vverlet: min_image_nopbc vverlet harmonic_force nlist util
	$(MAKE_EXE)
	$(FINISH)


lj_vverlet: SUB_CFLAGS += -DLJ
lj_vverlet: SUFFIX = ljv
lj_vverlet: PARAMS=lj
lj_vverlet: min_image lj_force vverlet nlist util
	$(MAKE_EXE)
	$(FINISH)

harmonic_vverlet_anderson: SUB_CFLAGS += -DHARMONIC -DNO_PBC -DANDERSON
harmonic_vverlet_anderson: SUFFIX = hva
harmonic_vverlet_anderson: PARAMS=harmonic
harmonic_vverlet_anderson: min_image_nopbc anderson_thermostat vverlet harmonic_force nlist util
	$(MAKE_EXE)
	$(FINISH)

lj_vverlet_bussi: SUB_CFLAGS += -DLJ -DBUSSI
lj_vverlet_bussi: SUFFIX = ljvb
lj_vverlet_bussi: PARAMS=lj
lj_vverlet_bussi: min_image bussi_thermostat lj_force vverlet nlist util
	$(MAKE_EXE)
	$(FINISH)


lj_vverlet_anderson: SUB_CFLAGS += -DLJ -DANDERSON
lj_vverlet_anderson: SUFFIX = ljvb
lj_vverlet_anderson: PARAMS=lj
lj_vverlet_anderson: min_image anderson_thermostat lj_force vverlet nlist util
	$(MAKE_EXE)
	$(FINISH)


harmonic_vverlet_bussi: SUB_CFLAGS += -DHARMONIC -DBUSSI -DNO_PBC
harmonic_vverlet_bussi: SUFFIX = hvb
harmonic_vverlet_bussi: PARAMS=harmonic
harmonic_vverlet_bussi: min_image_nopbc bussi_thermostat harmonic_force vverlet util nlist
	$(MAKE_EXE)
	$(FINISH)


clean: 
	rm -f *.o
	rm -f *.gcda

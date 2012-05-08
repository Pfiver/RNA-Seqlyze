# Makefile for homologyteams
#

COMP = gcc

FILES = \
        algorithm.o \
        filter.o \
        list.o \
        output.o \
        read.o \
        homologyteams.c

LIBS = 

CFLAGS = -c  -o $@


homologyteams: $(FILES)
	$(COMP) -o homologyteams $(FILES) $(LIBS)

# ---- alphabetized list of modules ----
algorithm.o: algorithm.c algorithm.h filter.h globaldata.h list.h types.h
	$(COMP) $(CFLAGS) algorithm.c

filter.o: filter.c algorithm.h filter.h types.h output.h
	$(COMP) $(CFLAGS) filter.c

list.o: list.c list.h
	$(COMP) $(CFLAGS) list.c

output.o: output.c output.h algorithm.h types.h
	$(COMP) $(CFLAGS) output.c

read.o: read.c globaldata.h types.h list.h read.h
	$(COMP) $(CFLAGS) read.c

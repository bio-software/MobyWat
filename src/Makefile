# MobyWat
# Calculation of hydration structures of molecular surfaces and interfaces
#
# Programmed by C. Hetenyi
# Web www.mobywat.com
# E-mail questions@mobywat.com
#
# Last modified 14-06-2015

CC = gcc
CFLAGS = -c -g -Wall 
LDFLAGS =
LIBFLAGS = -lm

SOURCES = analysis.c \
          analysis.h \
          calcs.c \
          cluster_by_olist.c \
          cluster_by_olist.h \
          consensus_plists.c \
          consensus_plists.h \
          editing.c \
          editing.h \
          global.h \
          input.c \
          input.h \
          insertion.c \
          macros.h \
          make_lists.c \
          make_pool.c \
          match.c \
          match_eval.c \
          merge_plists.c \
          mobywat.c \
          mobywat.h \
          netdraw.c \
          netdraw.h \
          pdb2xtc.c \
          prediction.c \
          prediction.h \
          print.c \
          print.h \
          process_input.c \
          sel_ifsf_wats.c \
          xdrfile.c \
          xdrfile.h \
          xdrfile_xtc.c \
          xdrfile_xtc.h

OBJECTS = $(SOURCES:.c=.o)

all: mobywat
	
mobywat: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< $(LIBFLAGS) -o $@

clean:
	rm *.o mobywat

install:
	cp mobywat $(HOME)/bin



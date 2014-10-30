####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = mpicc -icc
CFLAGS = -Wall -g -O1 $(PAPI_INC) 
LIBS = -lm $(PAPI_LIB)

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c vol2mesh.h
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

binconv: binconv.c
	gcc -std=gnu99 -Wall -g -O0 binconv.c -o binconv

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg binconv $(LIBPOS)

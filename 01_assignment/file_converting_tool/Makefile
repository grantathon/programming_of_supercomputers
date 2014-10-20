####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = gcc
CFLAGS = -Wall -g -O0 
LIBS =

LIBPOS =
AR =
ARFLAGS = 

SRCS = convert_text_to_binary.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: convert_text_to_binary 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o convert_text_to_binary $(LIBPOS)

CC = gcc
CFLAGS = -std=gnu99 -Wall -g -O0 
LIBS =

LIBPOS =
AR =
ARFLAGS = 

SRCS = binconv.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: binconv 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o binconv $(LIBPOS)

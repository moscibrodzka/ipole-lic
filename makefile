CC = gcc
CFLAGS =  -fopenmp -I/usr/include -O3
LDFLAGS = -lm 

SRCIPO = \
lic.c image.c

OBJIPO = \
lic.o image.o 

ipole: $(OBJIPO) makefile 
	$(CC) $(CFLAGS) -o lic $(OBJIPO) $(LDFLAGS)

$(OBJIPO) : makefile decs.h

clean:
	rm *.o *.ppm




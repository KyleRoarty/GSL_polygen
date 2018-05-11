CFLAGS= -g -Wall -std=gnu11
LDLIBS=
#`-Wl,-rpath -Wl,/usr/local/lib needed for gsl to work
LDFLAGS= -Wl,-rpath -Wl,/usr/local/lib -lgsl -lgslcblas
CC=gcc

SRC=./src
BIN=./bin

$(BIN)/$(P): $(SRC)/$(P).c
	$(CC) $(CFLAGS) $< $(LDLIBS) $(LDFLAGS) -o $@

clean:
	rm $(BIN)/*

CFLAGS= -g -Wall -std=gnu11
LDLIBS=
#`-Wl,-rpath -Wl,/usr/local/lib needed for gsl to work
LDFLAGS= -Wl,-rpath -Wl,/usr/local/lib -lgsl -lgslcblas -lm
CC=gcc

SRC=./src
BIN=./bin

$(BIN)/$(P): $(SRC)/$(P).c | $(BIN)
	$(CC) $(CFLAGS) $< $(LDLIBS) $(LDFLAGS) -o $@

$(BIN):
	mkdir -p $@

.PHONY: clean
clean:
	-rm -f $(BIN)/*


# Could also use 'pkg-config --cflags glib-2.0'
CFLAGS= -g -Wall -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -std=gnu11
LDLIBS=
# Could also use 'pkg-config --libs glib-2.0'
#`-Wl,-rpath -Wl,/usr/local/lib needed for gsl to work
LDFLAGS= -Wl,-rpath -Wl,/usr/local/lib -lgsl -lgslcblas -lm -lglib-2.0
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


INCLUDE_DIR = ./include
SRC_DIR = ./src
SHELL := /bin/bash
LIBS = -lgsl -lgslcblas -lm

CC = gcc
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall -Werror -Wpedantic
TARGET = TIDES

CCR = icc
CFLAGSR = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall -diag-disable=10441
TARGETR = TIDESR

DEPENDENCIES = $(SRC_DIR)/algelin3d.c $(SRC_DIR)/celmec.c $(SRC_DIR)/dynamical_system.c $(SRC_DIR)/parsing.c

.PHONY: example clean
.SILENT: tides tidesr clean

tides:
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

tidesr:
	source $$MODULESHOME/init/bash; \
	module --quiet load icc/latest; \
	module load gsl/2.7.1; \
	$(CCR) $(CFLAGSR) -o $(TARGETR) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

run:
	./$(TARGET) $(input)

plot:
	./$(TARGET) $(input) plot

example: tides
	./$(TARGET) example/EM_simplest_scenario.dat
	./$(TARGET) example/EM_simplest_scenario.dat plot

clean:
	-rm -f $(TARGET)

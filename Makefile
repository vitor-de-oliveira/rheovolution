INCLUDE_DIR = ./include
SRC_DIR = ./src
SHELL := /bin/bash
LIBS = -lgsl -lgslcblas -lm
TARGET = TIDES

CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall -diag-disable=10441
CC = icc

CFLAGSL = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall
CCL = gcc

DEPENDENCIES = $(SRC_DIR)/dynamical_system.c

.PHONY: cluster local clean
.SILENT: cluster local clean

cluster:
	source $$MODULESHOME/init/bash; \
	module --quiet load icc/latest; \
	module load gsl/2.7.1; \
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

local:
	$(CCL) $(CFLAGSL) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

clean:
	-rm -f $(TARGET)

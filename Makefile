INCLUDE_DIR = ./include
SRC_DIR = ./src
SHELL := /bin/bash
LIBS = -lgsl -lgslcblas -lm

TARGET = TIDES

CC = gcc
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall -Werror -Wpedantic

DEPENDENCIES = $(SRC_DIR)/algelin3d.c $(SRC_DIR)/celmec.c $(SRC_DIR)/dynamical_system.c $(SRC_DIR)/parsing.c

.PHONY: tides_remote example clean clean_slurm compile_icc
.SILENT: tides tides_remote clean clean_slurm compile_icc

tides_remote: compile_icc tides

tides:
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

run:
	./$(TARGET) $(input)

plot:
	./$(TARGET) $(input) plot

example: tides
	./$(TARGET) example/EM_simplest_scenario.dat
	./$(TARGET) example/EM_simplest_scenario.dat plot

clean:
	-rm -f $(TARGET)

clean_slurm:
	-rm -f slurm-*.out

compile_icc:
	$(eval CC = icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall -Werror -diag-disable=10441)
MAKEFLAGS += --no-print-directory
INCLUDE_DIR = ./include
SRC_DIR = ./src
SHELL := /bin/bash
LIBS = -lgsl -lgslcblas -lm

TARGET = TIDES

CC = cc
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall -Werror -Wpedantic

DEPENDENCIES = $(SRC_DIR)/algelin3d.c $(SRC_DIR)/celmec.c $(SRC_DIR)/dynamical_system.c $(SRC_DIR)/parsing.c

.PHONY: run orbital orientation plot all example_1 example_2 examples clean clean_slurm check_input
.SILENT: tides clean clean_slurm

tides:
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 -march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

run: check_input tides
	./$(TARGET) $(INPUT)

orbital: check_input tides
	./$(TARGET) $(INPUT) orbital

orientation: check_input tides
	./$(TARGET) $(INPUT) orientation

plot: check_input tides
	./$(TARGET) $(INPUT) plot

all: check_input tides run orbital orientation plot

example_1:
	$(eval EXAMPLE_NAME := Earth-Moon system)
	$(eval EXAMPLE_DIR := examples/EM_example.dat)
	@echo "Running $(EXAMPLE_NAME) example."
	@$(MAKE) all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME) example."

example_2:
	$(eval EXAMPLE_NAME := Earth-Moon-Sun system)
	$(eval EXAMPLE_DIR := examples/EMS_example.dat)
	@echo "Running $(EXAMPLE_NAME) example."
	@$(MAKE) all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME) example."

examples: example_1 example_2

clean:
	-rm -f $(TARGET)

clean_slurm:
	-rm -f slurm-*.out

check_input:
ifeq ($(INPUT), )
	@echo "Please, provide the input file via INPUT=<input_file_path>."
	@false
endif
MAKEFLAGS += --no-print-directory
INCLUDE_DIR = ./include
SRC_DIR = ./src
SHELL := /bin/bash
LIBS = -lgsl -lgslcblas -lm

TARGET = TIDES

CC = cc
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
		 -march=native -Wall -Werror -Wpedantic \
#		 -fsanitize=address

DEPENDENCIES =  $(SRC_DIR)/linear_algebra.c \
				$(SRC_DIR)/celestial_mechanics.c \
				$(SRC_DIR)/tidal_theory.c \
				$(SRC_DIR)/dynamical_system.c \
				$(SRC_DIR)/data_processing.c

.PHONY: run orbit spin plot all example_1 example_2 example_3 examples \
		clean clean_slurm input_check
.SILENT: tides calibration clean clean_slurm

tides:
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

run: input_check tides
	./$(TARGET) $(INPUT)

orbit: input_check tides
	./$(TARGET) $(INPUT) orbit

spin: input_check tides
	./$(TARGET) $(INPUT) spin

plot: input_check tides
	./$(TARGET) $(INPUT) plot

all: input_check tides run orbit spin plot

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

example_3:
	$(eval EXAMPLE_NAME := rigid Earth)
	$(eval EXAMPLE_DIR := examples/E_example.dat)
	@echo "Running $(EXAMPLE_NAME) example."
	@$(MAKE) all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME) example."

examples: example_1 example_2 example_3

calibration:
	$(eval TARGET = CALIBRATION)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
	./$(TARGET)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/calibration.c $(DEPENDENCIES) $(LIBS)
	./$(TARGET)

clean:
	-rm -f $(TARGET)

clean_slurm:
	-rm -f slurm-*.out

input_check:
ifeq ($(INPUT), )
	@echo "Please, provide the input file via INPUT=<input_file_path>."
	@false
endif
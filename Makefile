MAKEFLAGS += --no-print-directory
INCLUDE_DIR = ./include
SRC_DIR = ./src
SHELL := /bin/bash
LIBS = -lgsl -lgslcblas -lm

TARGET = tides

CC = cc
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
		 -march=native -Wall -Werror -Wpedantic \
#		 -fsanitize=address

DEPENDENCIES =  $(SRC_DIR)/linear_algebra.c \
				$(SRC_DIR)/celestial_mechanics.c \
				$(SRC_DIR)/tidal_theory.c \
				$(SRC_DIR)/dynamical_system.c \
				$(SRC_DIR)/data_processing.c

.PHONY: compile run orbit spin plot all example_1 example_2 example_3 \
    examples examples_parallel tests tests_parallel run_calibration \
    clean clean_calibration clean_examples clean_slurm input_check
.SILENT: tides compile examples_parallel tests_parallel clean clean_examples clean_slurm

tides: $(SRC_DIR)/main.c $(DEPENDENCIES)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

compile: $(SRC_DIR)/main.c $(DEPENDENCIES)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/main.c $(DEPENDENCIES) $(LIBS)

run: input_check
	./$(TARGET) $(INPUT)

orbit: input_check
	./$(TARGET) $(INPUT) orbit

spin: input_check
	./$(TARGET) $(INPUT) spin

plot: input_check
	./$(TARGET) $(INPUT) plot

all: input_check run orbit spin plot

example_E:
	$(eval EXAMPLE_NAME := Earth)
	$(eval EXAMPLE_DIR := examples/E_example.dat)
	@echo "Running $(EXAMPLE_NAME) example."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME) example."

example_EM:
	$(eval EXAMPLE_NAME := Earth-Moon system)
	$(eval EXAMPLE_DIR := examples/EM_example.dat)
	@echo "Running $(EXAMPLE_NAME) example."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME) example."

example_EMS:
	$(eval EXAMPLE_NAME := Earth-Moon-Sun system)
	$(eval EXAMPLE_DIR := examples/EMS_example.dat)
	@echo "Running $(EXAMPLE_NAME) example."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME) example."

examples: example_E example_EM example_EMS

examples_parallel:
	@$(MAKE) -j examples

test_1:
	$(eval EXAMPLE_NAME := test 1: rigid Earth)
	$(eval EXAMPLE_DIR := examples/test_1.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

test_2:
	$(eval EXAMPLE_NAME := test 2: Earth-Moon system)
	$(eval EXAMPLE_DIR := examples/test_2.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

test_3:
	$(eval EXAMPLE_NAME := test 3: Earth-Moon-Sun system)
	$(eval EXAMPLE_DIR := examples/test_3.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

tests: test_1 test_2 test_3

tests_parallel:
	@$(MAKE) -j tests

calibration: $(SRC_DIR)/calibration.c $(DEPENDENCIES)
	$(eval TARGET = calibration)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
	./$(TARGET)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/calibration.c $(DEPENDENCIES) $(LIBS)

run_calibration:
	$(eval TARGET = calibration)
	./$(TARGET)

clean:
	-rm -f $(TARGET)

clean_calibration:
	$(eval TARGET = calibration)
	-rm -f $(TARGET)

clean_examples:
	-rm -f -r examples/output_*/

clean_slurm:
	-rm -f slurm-*.out

input_check:
ifeq ($(INPUT), )
	@echo "Please, provide the input file via INPUT=<input_file_path>."
	@false
endif
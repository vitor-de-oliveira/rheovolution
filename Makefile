MAKEFLAGS += --no-print-directory
INCLUDE_DIR = ./include
SRC_DIR = ./src
SHELL := /bin/bash
LIBS = -l:libgsl.so.27 -l:libgslcblas.so.0 -lm

TARGET = rheo

MAIN = main.c

CC = cc
CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
		 -march=native -Wall -Werror -Wpedantic

DEPENDENCIES =  $(SRC_DIR)/linear_algebra.c \
				$(SRC_DIR)/celestial_mechanics.c \
				$(SRC_DIR)/tidal_theory.c \
				$(SRC_DIR)/dynamical_system.c \
				$(SRC_DIR)/data_processing.c

.PHONY: compile run orbit spin plot all example_1_EMS \ 
	example_2_E_rigid example_3_E_deformable example_4_EM_Maxwell \ 
	example_5_EM_Burgers_wobble example_6_EM_Burgers_drift \
    examples examples_parallel tests tests_parallel \
    clean clean_examples clean_slurm input_check
.SILENT: rheo compile examples_parallel tests_parallel clean \ 
	clean_examples clean_slurm

rheo: $(SRC_DIR)/$(MAIN) $(DEPENDENCIES)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/$(MAIN) $(DEPENDENCIES) $(LIBS)

compile: $(SRC_DIR)/$(MAIN) $(DEPENDENCIES)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC_DIR)/$(MAIN) $(DEPENDENCIES) $(LIBS)

debug: $(SRC_DIR)$(MAIN) $(DEPENDENCIES)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -g -o $(TARGET) $(SRC_DIR)/$(MAIN) $(DEPENDENCIES) $(LIBS)

sanitize: $(SRC_DIR)/$(MAIN) $(DEPENDENCIES)
ifeq ($(CC),icc)
	$(eval CFLAGS = -std=c11 -I$(INCLUDE_DIR) -D_XOPEN_SOURCE -O3 \
					-march=native -Wall -Werror -diag-disable=10441)
endif
	$(CC) $(CFLAGS) -fsanitize=address -o $(TARGET) $(SRC_DIR)/$(MAIN) $(DEPENDENCIES) $(LIBS)

run: input_check
	./$(TARGET) $(INPUT)

orbit: input_check
	./$(TARGET) $(INPUT) orbit

spin: input_check
	./$(TARGET) $(INPUT) spin

plot: input_check
	./$(TARGET) $(INPUT) plot

all: input_check run orbit spin plot

example_1_EMS:
	$(eval EXAMPLE_NAME := Example 1: Earth-Moon-Sun system)
	$(eval EXAMPLE_DIR := examples/Example_1_EMS.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

example_2_E_rigid:
	$(eval EXAMPLE_NAME := Example 2: Rigid Earth)
	$(eval EXAMPLE_DIR := examples/Example_2_E_rigid.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

example_3_E_deformable:
	$(eval EXAMPLE_NAME := Example 3: Deformable Earth with Maxwell rheology - Chandler Wobble)
	$(eval EXAMPLE_DIR := examples/Example_3_E_deformable.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

example_4_EM_Maxwell:
	$(eval EXAMPLE_NAME := Example 4: Earth-Moon system with Maxwell rheology - lunar drift)
	$(eval EXAMPLE_DIR := examples/Example_4_EM.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

example_5_EM_Burgers_wobble:
	$(eval EXAMPLE_NAME := Example 5: Earth-Moon system with Burgers rheology - Chandler Wobble)
	$(eval EXAMPLE_DIR := examples/Example_5_EM_Burgers_Chandler_Wobble.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

example_6_EM_Burgers_drift:
	$(eval EXAMPLE_NAME := Example 6: Earth-Moon system with Burgers rheology - lunar drift)
	$(eval EXAMPLE_DIR := examples/Example_6_EM_Burgers_Moon_orbital_drift.dat)
	@echo "Running $(EXAMPLE_NAME)."
	@$(MAKE) -j 1 all INPUT=$(EXAMPLE_DIR)
	@echo "Finished running $(EXAMPLE_NAME)."

examples: example_1_EMS example_2_E_rigid example_3_E_deformable example_4_EM_Maxwell example_5_EM_Burgers_wobble example_6_EM_Burgers_drift

examples_parallel:
	@$(MAKE) -j examples

clean:
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

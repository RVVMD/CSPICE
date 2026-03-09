# MNA Solver Makefile
# Simple alternative to CMake

CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -fPIC
LDFLAGS = -lm

# Directories
SRC_DIR = src
ELEM_DIR = elements
NONLIN_DIR = elements/nonlinear
INC_DIR = include

# Source files
SRCS = $(SRC_DIR)/matrix.c \
       $(SRC_DIR)/solver/core.c \
       $(SRC_DIR)/solver/dc.c \
       $(SRC_DIR)/solver/ac.c \
       $(SRC_DIR)/solver/transient.c \
       $(ELEM_DIR)/passive.c \
       $(ELEM_DIR)/sources.c \
       $(NONLIN_DIR)/nonlinear.c

# Object files
OBJS = $(SRCS:.c=.o)

# Include paths
INCLUDES = -I. -I$(INC_DIR) -I$(ELEM_DIR) -I$(SRC_DIR)

# Targets
LIB_STATIC = libmna_solver.a
LIB_SHARED = libmna_solver.so

.PHONY: all clean static shared

all: static shared

# Static library
static: $(LIB_STATIC)

$(LIB_STATIC): $(OBJS)
	ar rcs $@ $^

# Shared library
shared: $(LIB_SHARED)

$(LIB_SHARED): $(OBJS)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

# Compile source files
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Example executable (if TR5_sim.c exists)
tr5_sim: TR5_sim.c $(LIB_STATIC)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIB_STATIC) $(LDFLAGS)

clean:
	rm -f $(OBJS) $(LIB_STATIC) $(LIB_SHARED) tr5_sim

# Dependencies (header files)
$(SRC_DIR)/matrix.o: $(INC_DIR)/types.h $(INC_DIR)/matrix.h
$(SRC_DIR)/solver/core.o: $(INC_DIR)/types.h $(INC_DIR)/matrix.h $(SRC_DIR)/solver/core.h
$(SRC_DIR)/solver/dc.o: $(INC_DIR)/types.h $(INC_DIR)/matrix.h $(SRC_DIR)/solver/dc.h \
                        $(ELEM_DIR)/passive.h $(ELEM_DIR)/sources.h $(NONLIN_DIR)/nonlinear.h
$(SRC_DIR)/solver/ac.o: $(INC_DIR)/types.h $(INC_DIR)/matrix.h $(SRC_DIR)/solver/ac.h \
                        $(NONLIN_DIR)/nonlinear.h
$(SRC_DIR)/solver/transient.o: $(INC_DIR)/types.h $(INC_DIR)/matrix.h $(SRC_DIR)/solver/transient.h \
                               $(ELEM_DIR)/passive.h $(ELEM_DIR)/sources.h $(NONLIN_DIR)/nonlinear.h
$(ELEM_DIR)/passive.o: $(INC_DIR)/types.h $(ELEM_DIR)/passive.h $(INC_DIR)/matrix.h
$(ELEM_DIR)/sources.o: $(INC_DIR)/types.h $(ELEM_DIR)/sources.h $(INC_DIR)/matrix.h
$(NONLIN_DIR)/nonlinear.o: $(INC_DIR)/types.h $(NONLIN_DIR)/nonlinear.h $(INC_DIR)/matrix.h

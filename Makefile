# compiler
FC = gfortran
CC = gcc

# compile flags
FCFLAGS = -g -Ofast -ffixed-line-length-0 -std=legacy -Iinclude
CFLAGS = -g -Ofast -Wno-pointer-to-int-cast -Iinclude

SRC_DIR = src
OBJ_DIR = lib
BIN_DIR = .

# source files and objects
FSRCS = $(wildcard $(SRC_DIR)/*.f90)
FFSRCS = $(wildcard $(SRC_DIR)/*.F)
CSRCS = $(wildcard $(SRC_DIR)/*.c)

FOBJ = $(FSRCS:$(SRC_DIR)/%.f90=$(OBJ_DIR)/%.o)
FFOBJ = $(FFSRCS:$(SRC_DIR)/%.F=$(OBJ_DIR)/%.o)
COBJ = $(CSRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# program name
PROGRAM = onepigen

MKDIR_P = mkdir -p
.PHONY: directories


all: directories spp_tbl $(COBJ) $(FOBJ) $(FFOBJ) $(PROGRAM)

spp_tbl:
	tar -xvf spp_tbl.tar.gz

directories: ${OBJ_DIR} ${BIN_DIR}

${OBJ_DIR}:
	${MKDIR_P} ${OBJ_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

$(PROGRAM): $(FOBJ) $(COBJ)
	$(FC) $(FCFLAGS) -o ${BIN_DIR}/$@_lund $(COBJ) $(FOBJ) $(FFOBJ)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F
	$(FC) $(FCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	rm -f ${OBJ_DIR}/*.o ${BIN_DIR}/$(PROGRAM)_lund

##
## @file		makefile
## @brief		PICSP makefile.
## @author		Sayan Adhikari <sayan.adhikari@fys.uio.no>
##

CC		= g++
COPT	= -g -O3
DOPT 	= -O0

CLOCAL = 	-Ilib/iniparser/src\
			-lm
CLOCALSA = 	-I /opt/local/include\
			
LLOCAL =	-Ilib/iniparser/src\
			-lm
LLOCALSA =      -I /opt/local/include\
                        
-include local.mk

EXEC	= picsp
CADD	= # Additional CFLAGS accessible from CLI
CFLAGS	= -std=c11 -Wall $(CLOCAL) $(CLOCALSA) $(COPT) $(CADD) # Flags for compiling

LFLAGS	= -std=c11 -Wall $(LLOCAL) $(LLOCALSA) $(COPT) $(CADD) # Flags for linking

SDIR	= src
ODIR	= src/obj
HDIR	= src
LDIR	= lib


HEAD_	=
SRC_	=

OBJ_	= $(SRC_:.c=.o)


HEAD	= $(patsubst %,$(HDIR)/%,$(HEAD_))
SRC		= $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ		= $(patsubst %,$(ODIR)/%,$(OBJ_))




LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h

LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

all: version $(EXEC)

local: version $(EXEC).local


$(EXEC).local: $(ODIR)/main.local.o $(OBJ) $(LIBOBJ)
	@echo "Linking PICSP (using main.local.c)"
	@$(CC) $^ -o $(EXEC) $(LFLAGS)
	@echo "PICSP is built"

$(EXEC): $(ODIR)/main.o $(OBJ) $(LIBOBJ)
	@echo "Linking PICSP"
	@$(CC) $^ -o $@ $(LFLAGS)
	@echo "PICSP is built"

$(ODIR)/%.o: $(SDIR)/%.c $(HEAD)
	@echo "Compiling $<"
	@echo $(HEAD) | xargs -n1 ./aux/check.sh
	@mkdir -p $(ODIR)
	@./aux/check.sh $<
	@$(CC) -c $< -o $@ $(CFLAGS)


$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	@echo "Building iniparser"
	@cd $(LDIR)/iniparser && $(MAKE) libiniparser.a > /dev/null 2>&1
	
.phony: version
version:
	@echo "Embedding git version"
	@echo "#define VERSION \"$(shell git describe --abbrev=4 --dirty --always --tags)\"" > $(SDIR)/version.h

veryclean: clean
	@echo "Cleaning executable and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) veryclean > /dev/null 2>&1


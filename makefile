##
## @file		makefile
## @brief		PICSP makefile.
## @author		Sayan Adhikari <sayan.adhikari@fys.uio.no>
##

CC		= g++

ARG = -std=c++11 -Wall
                        

EXEC	= picsp

CFLAGS	=  -Ilib/iniparser/src   # Flags for compiling

LFLAGS	= -Llib/iniparser -liniparser # Flags for linking

SDIR	= src
ODIR	= src/obj
#HDIR	= src
LDIR	= lib


LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h
OBJ_	= $(SRC_:.c=.o)

OBJ		= $(patsubst %,$(ODIR)/%,$(OBJ_))
LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

all: version $(EXEC)


$(ODIR)/%.o: $(SDIR)/%.c
	@echo "Compiling $<"
	@$(CC) $(ARG) -c $< -o $@ $(CFLAGS)

$(EXEC): $(OBJ) $(LIBOBJ)
	@echo "Linking PICSP"
	@$(CC) $(ARG) $^ -o $@ $(LFLAGS)
	@echo "PICSP is built"

$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	@echo "Building iniparser"
	@cd $(LDIR)/iniparser && $(MAKE) libiniparser.a > /dev/null 2>&1
	
.phony: version
version:
	@echo "Embedding git version"
	@echo "#define VERSION \"$(shell git describe --abbrev=4 --dirty --always --tags)\"" > $(SDIR)/version.h

clean:
	@echo "Cleaning compilation files (run \"make veryclean\" to clean more)"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o

veryclean: clean
	@echo "Cleaning executable and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) veryclean > /dev/null 2>&1


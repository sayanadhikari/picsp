##
## @file		makefile
## @brief		PICSP makefile.
## @author		Sayan Adhikari <sayan.adhikari@fys.uio.no>
##

CXX		= g++

#ARG =

EXEC	= picsp

CXXFLAGS = -g -std=c++11 -Wall
CFLAGS	=  -Ilib/iniparser/src   # Flags for compiling

LFLAGS	=  -Llib/iniparser -liniparser # Flags for linking

SDIR	= src
ODIR	= src/obj
#HDIR	= src
LDIR	= lib
OUTDIR  = output

LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h

LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

all: version $(EXEC)

OBJ = $(ODIR)/main.o
SRC = $(SDIR)/main.cpp

$(EXEC): $(OBJ)
	@echo "PICSP is being compiled"
	@mkdir -p $(ODIR)
	@$(CXX) $(CXXFLAGS) -o  $(EXEC) $(OBJ) $(CFLAGS) $(LFLAGS)
	@echo "PICSP is built"

$(OBJ): $(SRC)
	@$(CXX) $(CXXFLAGS) -c $(SRC) -o $(OBJ) $(CFLAGS) $(LFLAGS)

$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	@echo "Building iniparser"
	@cd $(LDIR)/iniparser && $(MAKE) libiniparser.a > /dev/null 2>&1

.phony: version
version:
	@echo "Embedding git version"
	@echo "#define VERSION \"$(shell git describe --abbrev=4 --dirty --always --tags)\"" > $(SDIR)/version.h

clean:
	@echo "Cleaning compiled files"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o
	@rm -f *.dat
	@rm -rf $(OUTDIR)

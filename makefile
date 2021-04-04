##
## @file		makefile
## @brief		PICSP makefile.
## @author		Sayan Adhikari <sayan.adhikari@fys.uio.no>
##

CXX		= g++

CXXLOCAL = -Ilib/iniparser/src
LLOCAL = -Ilib/iniparser/src

FFLAGS = -lfftw3 -lm
HFLAGS = -lhdf5 -lhdf5_cpp

EXEC	= picsp

CXXFLAGS = -g -std=c++11 -Wall $(CXXLOCAL)  # Flags for compiling
LFLAGS	=  -g -std=c++11 -Wall $(LLOCAL) # Flags for linking

SDIR	= src
ODIR	= src/obj
#HDIR	= src
LDIR	= lib
OUTDIR  = output

SRC_ 	= # Additional CPP files
OBJ_	= $(SRC_:.cpp=.o)

SRC = $(patsubst %,$(SDIR)/%,$(SRC_))
OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))

LIBOBJ_	= iniparser/libiniparser.a
LIBHEAD_= iniparser/src/iniparser.h


LIBOBJ = $(patsubst %,$(LDIR)/%,$(LIBOBJ_))
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))


all: version $(EXEC)

$(EXEC): $(ODIR)/main.o $(OBJ) $(LIBOBJ)
	@echo "Linking PICSP"
	@$(CXX) $^ -o $@ $(LFLAGS) $(FFLAGS) $(HFLAGS)
	@echo "PICSP is built"

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $<"
	@mkdir -p $(ODIR)
	@mkdir -p $(OUTDIR)
	@$(CXX) -c $< -o $@ $(CXXFLAGS)

$(LDIR)/iniparser/libiniparser.a: $(LIBHEAD)
	@echo "Building iniparser"
	@cd $(LDIR)/iniparser && $(MAKE) libiniparser.a > /dev/null 2>&1

.phony: version
version:
	@echo "Embedding git version"
	@echo "#define VERSION \"$(shell git describe --abbrev=4 --dirty --always --tags)\"" > $(SDIR)/version.h

clean:
	@echo "Cleaning compiled files. (run 'make veryclean' to remove executables and more)"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o $(SDIR)/*~
	@rm -f *.dat
	@rm -rf $(OUTDIR)
veryclean: clean
	@echo "Cleaning executables and iniparser"
	@rm -f $(EXEC)
	@cd $(LDIR)/iniparser && $(MAKE) veryclean > /dev/null 2>&1

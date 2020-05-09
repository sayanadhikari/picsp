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

OBJ     = picsp.o

all: iniparse version $(EXEC)

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) $(ARG) -c -o $@ $< $(CFLAGS) $(LFLAGS)

$(EXEC): $(OBJ)
	@echo "Compiling & Linking PICSP"
	@$(CC) $(ARG) -o $@ $^ $(CFLAGS) $(LFLAGS)
	@echo "PICSP is built"

iniparse:
	@echo "Building iniparser"
    @cd $(LDIR)/iniparser && @make

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


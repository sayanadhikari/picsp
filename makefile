##
## @file		makefile
## @brief		PICSP makefile.
## @author		Sayan Adhikari <sayan.adhikari@fys.uio.no>
##

CC		= g++ -std=c++11 -Wall

#ARG =
                        

EXEC	= picsp

CFLAGS	=  -Ilib/iniparser/src   # Flags for compiling

LFLAGS	=  -Llib/iniparser -liniparser # Flags for linking

SDIR	= src
ODIR	= src/obj
#HDIR	= src
LDIR	= lib

OBJ     = picsp.o

LIBHEAD_= iniparser/src/iniparser.h
LIBHEAD = $(patsubst %,$(LDIR)/%,$(LIBHEAD_))

all: version $(EXEC)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(LFLAGS)

$(EXEC): $(OBJ)
	@echo "Compiling & Linking PICSP"
	@$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)
	@echo "PICSP is built"


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


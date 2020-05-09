##
## @file		makefile
## @brief		PICSP makefile.
## @author		Sayan Adhikari <sayan.adhikari@fys.uio.no>
##

CC		= g++ -std=c++11 -Wall

#ARG =
                        

#EXEC	= picsp

CFLAGS	=  -Ilib/iniparser/src   # Flags for compiling

LFLAGS	=  -Llib/iniparser -liniparser # Flags for linking

SDIR	= src
ODIR	= src/obj
#HDIR	= src
LDIR	= lib
OUTDIR  = output

OBJ = picsp

all: version $(EXEC)

$(ODIR)/%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) $(LFLAGS)

$(EXEC): $(OBJ)
	@echo "Compiling & Linking PICSP"
	@$(CC) -c $(SDIR)/%.cpp -o $(OBJ) $(CFLAGS) $(LFLAGS)
	@echo "PICSP is built"


.phony: version
version:
	@echo "Embedding git version"
	@echo "#define VERSION \"$(shell git describe --abbrev=4 --dirty --always --tags)\"" > $(SDIR)/version.h

clean:
	@echo "Cleaning compilation files (run \"make veryclean\" to clean more)"
	@rm -f *~ $(ODIR)/*.o $(SDIR)/*.o
	@rm -f $(OUTDIR)/*



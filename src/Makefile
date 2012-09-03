# no console output during compilation:
#.SILENT:

# *******************************************************
# ***          Comecar por limpar os sufixos          ***
# *******************************************************
.SUFFIXES:

# *******************************************************
# ***     Especificar os sufixos para .f .o .do       ***
# *******************************************************
#.SUFFIXES: .f90 .o

#
# *******************************************************
# ***                       Macros                    ***
# *******************************************************
#FC = /usr/bin/f77
#FC	= g77 -Wall
FC	= gfortran
FFLAGS = -O0 -ffree-form -fno-align-commons
#dbx     = -O5 -r8 -g
#profil  = -p -O5 -r8 
#samedir = .
#FTN     = ftnchek

#
# *******************************************************
# *** Regra que por defeito produz os ficheiros .o **
# *******************************************************
%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $*.f
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $*.f90

#
# *******************************************************
# ***   Especificar as directorias com as subrotinas  ***
# *******************************************************
FNOSC = /home/joao/Documents/Codes/FNOSC/src
LIB = /home/joao/Documents/Codes/FNOSC/lib

# plot library:
#JFp = /home/joao/Utils/gnuFor2
# spline library:
#JFspl = /home/joao/Programs/CODE_freqFit/dierckx

#----------------------------------------------------------

FILES = \
$(FNOSC)/consts.o \
$(FNOSC)/commonvar.o \
$(FNOSC)/consl.o \
$(FNOSC)/fehl.o \
$(FNOSC)/goout0.o \
$(FNOSC)/gooutl.o \
$(FNOSC)/modeid.o \
$(FNOSC)/puls.o \
$(FNOSC)/quad.o \
$(FNOSC)/read.o \
$(FNOSC)/rkf.o \
$(FNOSC)/rkfcow.o \
$(FNOSC)/rkfrad.o \
$(FNOSC)/rkfs.o \
$(FNOSC)/splntl.o \
$(FNOSC)/read_sun.o \
$(FNOSC)/output.o

LIBS = \
$(LIB)/lib_array.o


#FPLOT = \
#$(JFp)/gnufor2.o

#FSPLINE = \
#$(JFspl)/concur.o \
#$(JFspl)/curev.o \
#$(JFspl)/fpadpo.o \
#$(JFspl)/fpback.o \
#$(JFspl)/fpbspl.o \
#$(JFspl)/fpched.o \
#$(JFspl)/fpcons.o \
#$(JFspl)/fpdisc.o \
#$(JFspl)/fpgivs.o \
#$(JFspl)/fpinst.o \
#$(JFspl)/fpknot.o \
#$(JFspl)/fppocu.o \
#$(JFspl)/fprati.o \
#$(JFspl)/fprota.o 


# **********************************************************
# ***             Compilar os programas                  *** 
# **********************************************************

puls: $(LIBS) $(FILES) 
	$(FC) $(FILES) $(LIBS) -o ../puls

clean:
	rm -f $(FILES) $(LIBS) *~ ../puls

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------


#CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Enterprise_Server1263/CPLEX_Sudio/cplex
#CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Enterprise_Server1263/CPLEX_Sudio/concert

CPLEXDIR      = /modules/apps/cplex/12.7/cplex
CONCERTDIR    = /modules/apps/cplex/12.7/concert

# ---------------------------------------------------------------------
# Compiler selection
# ---------------------------------------------------------------------

CCC = g++ -g
CC  = gcc -g

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------

COPT = -m64 -O3 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CPLEX_FLAGS = -lilocplex -lconcert -lcplex -lm -lpthread -fexceptions -lrt

.PHONY: all clean

CPLEXINCDIR   = $(CPLEXDIR)/include
CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)

RAIZ = /lustre/desymcb/PIDLS-NT/Raiz
BASE = /lustre/desymcb/PIDLS-NT/SCF_FL/Base
#------------------------------------------------------------
#  make all      : to compile the examples.
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX =  $(RAIZ)/raiz $(BASE)/scf

all: $(C_EX) main
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o scf_fl.o main.o -o scf $(CPLEX_FLAGS)
	clear
	qsub script

# ------------------------------------------------------------
# The examples
#

$(RAIZ)/raiz: $(RAIZ)/raiz.h $(RAIZ)/raiz.cpp
	$(CCC) -c $(CFLAGS) $(RAIZ)/raiz.cpp -o raiz.o $(CPLEX_FLAGS)

$(BASE)/scf: $(BASE)/scf_fl.h $(BASE)/scf_fl.cpp
	$(CCC) -c $(CFLAGS) $(BASE)/scf_fl.cpp -o scf_fl.o $(CPLEX_FLAGS)

main: main.cpp
	$(CCC) -c $(CFLAGS) main.cpp -o main.o $(CPLEX_FLAGS)
clean :
	rm -r *.o

# Local Variables:
# mode: makefile
# End:

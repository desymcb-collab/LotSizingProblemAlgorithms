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

CCOPT = -m64 -O3 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CPLEX_FLAGS = -lilocplex -lconcert -lcplex -lm -lpthread -fexceptions

.PHONY: all clean

CPLEXINCDIR   = $(CPLEXDIR)/include
CFLAGS  = $(CCOPT)  -I$(CPLEXINCDIR)


#SBG = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/Subgradiente
RAIZ = /home/desymcb/Desiree/PIDLS-NT/Raiz
#MO = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_SeLSP/NaoTriangular/MaquinasParalelas-SetupTri/SCF_FL
#MSG = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/ModelosSubGrad
#NF = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/NonaFactibiliza
#MF = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/MamisFactibiliza

#------------------------------------------------------------
#  make all      : to compile the examples.
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX =  $(RAIZ)/raiz scf_fl

all: $(C_EX) main
	$(CCC) $(CFLAGS) $(CCLNDIRS) $(RAIZ)/raiz.o scf_fl.o main.o -o scf_fl_exec7 $(CPLEX_FLAGS)
	clear
	qsub script

# ------------------------------------------------------------
# The examples
#

$(RAIZ)/raiz: $(RAIZ)/raiz.h $(RAIZ)/raiz.cpp
	$(CCC) -c $(CFLAGS) $(RAIZ)/raiz.cpp -o $(RAIZ)/raiz.o $(CPLEX_FLAGS)

scf_fl: scf_fl.h scf_fl.cpp
	$(CCC) -c $(CFLAGS) scf_fl.cpp -o scf_fl.o $(CPLEX_FLAGS)

main: main.cpp
	$(CCC) -c $(CFLAGS) main.cpp -o main.o $(CPLEX_FLAGS)
clean :
	rm -r *.o

# Local Variables:
# mode: makefile
# End:

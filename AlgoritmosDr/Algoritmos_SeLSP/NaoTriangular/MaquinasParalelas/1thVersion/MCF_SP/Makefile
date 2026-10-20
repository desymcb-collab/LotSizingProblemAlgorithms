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

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio126/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio126/concert

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
CPLEX_FLAGS = -lilocplex -lconcert -lcplex -lm -lpthread -fexceptions -lrt

.PHONY: all clean

CPLEXINCDIR   = $(CPLEXDIR)/include
CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)

#SBG = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/Subgradiente
RAIZ = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_SeLSP/NaoTriangular/MaquinasParalelas/Raiz
#MO = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_SeLSP/NaoTriangular/MaquinasParalelas-SetupTri/mcf_sp
#MSG = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/ModelosSubGrad
#NF = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/NonaFactibiliza
#MF = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_S&LSP/MamisFactibiliza

#------------------------------------------------------------
#  make all      : to compile the examples.
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX =  raiz mcf_sp

all: $(C_EX) main
	$(CCC) $(RAIZ)/raiz.o mcf_sp.o main.o -o mcf_sp_exec $(CPLEX_FLAGS)
	clear
	./mcf_sp_exec

# ------------------------------------------------------------
# The examples
#

raiz: $(RAIZ)/raiz.h $(RAIZ)/raiz.cpp
	$(CCC) -c $(CCOPT) $(RAIZ)/raiz.cpp -o $(RAIZ)/raiz.o $(CPLEX_FLAGS)

mcf_sp: mcf_sp.h mcf_sp.cpp
	$(CCC) -c $(CCOPT) mcf_sp.cpp -o mcf_sp.o $(CPLEX_FLAGS)

main: main.cpp
	$(CCC) -c $(CCOPT) main.cpp -o main.o $(CPLEX_FLAGS)
clean :
	rm -r *.o

# Local Variables:
# mode: makefile
# End:

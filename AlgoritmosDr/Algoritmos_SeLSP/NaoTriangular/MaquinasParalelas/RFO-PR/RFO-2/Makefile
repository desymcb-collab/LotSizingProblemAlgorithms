SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------


CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio126/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio126/concert

#CPLEXDIR      = /modules/apps/cplex/12.7/cplex
#CONCERTDIR    = /modules/apps/cplex/12.7/concert

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
CPLEX_FLAGS = -lilocplex -lconcert -lcplex -lm -lpthread -fexceptions -ldl

.PHONY: all clean

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CFLAGS  = $(COPT)  -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

RAIZ = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_SeLSP/NaoTriangular/MaquinasParalelas/Raiz
BASE = /home/desiree/Dropbox/Doutorado/AlgoritmosDr/Algoritmos_SeLSP/NaoTriangular/MaquinasParalelas/RFO-PR/Base
#------------------------------------------------------------
#  make all      : to compile the examples.
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX =  $(RAIZ)/raiz $(BASE)/rf

all: $(C_EX) main
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o main.o -o rfo $(CPLEX_FLAGS)
	clear
	./rfo
	rm -r *.o

# ------------------------------------------------------------
# The examples
#

$(RAIZ)/raiz: $(RAIZ)/raiz.h $(RAIZ)/raiz.cpp
	$(CCC) -c $(CFLAGS) $(CCLNDIRS) $(RAIZ)/raiz.cpp -o raiz.o $(CPLEX_FLAGS)

$(BASE)/rf: $(BASE)/rf.h $(BASE)/rf.cpp
	$(CCC) -c $(CFLAGS) $(CCLNDIRS) $(BASE)/rf.cpp -o rf.o $(CPLEX_FLAGS)

main: main.cpp
	$(CCC) -c $(CFLAGS) $(CCLNDIRS) main.cpp -o main.o $(CPLEX_FLAGS)
clean :
	rm -r *.o

# Local Variables:
# mode: makefile
# End:

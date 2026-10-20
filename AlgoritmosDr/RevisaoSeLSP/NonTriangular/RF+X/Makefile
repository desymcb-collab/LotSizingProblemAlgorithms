SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------


CPLEXDIR      = /home/desymcb/ibm/ILOG/CPLEX_Studio1210/cplex
CONCERTDIR    = /home/desymcb/ibm/ILOG/CPLEX_Studio1210/concert

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
CPLEX_FLAGS = -lilocplex -lconcert -lcplex -lm -lpthread -fexceptions -lrt -ldl

.PHONY: all clean

CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR   = $(CONCERTDIR)/include
CFLAGS  = $(COPT)  -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

#------------------------------------------------------------
#  make all      : to compile the examples.
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX =  raiz rf fo ks rfks

all: $(C_EX) main
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o fo.o ks.o rfks.o main.o -o HeurX $(CPLEX_FLAGS)
	rm -r *.o
	clear
	qsub script1
	qsub script2
	qsub script3
	qsub script4
	qsub script5
	qsub script6
	qsub script7
	qsub script8
	qsub script9
	qsub script10
	

# ------------------------------------------------------------
# The examples
#

raiz: raiz.h raiz.cpp
	$(CCC) -c $(CFLAGS) raiz.cpp -o raiz.o $(CPLEX_FLAGS)

rf: rf.h rf.cpp
	$(CCC) -c $(CFLAGS) rf.cpp -o rf.o $(CPLEX_FLAGS)

fo: fo.h fo.cpp
	$(CCC) -c $(CFLAGS) fo.cpp -o fo.o $(CPLEX_FLAGS)

ks: ks.h ks.cpp
	$(CCC) -c $(CFLAGS) ks.cpp -o ks.o $(CPLEX_FLAGS)

rfks: rfks.h rfks.cpp
	$(CCC) -c $(CFLAGS) rfks.cpp -o rfks.o $(CPLEX_FLAGS)

main: main.cpp
	$(CCC) -c $(CFLAGS) main.cpp -o main.o $(CPLEX_FLAGS)
clean :
	rm -r *.o

# Local Variables:
# mode: makefile
# End:

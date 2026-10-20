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

#CPLEXDIR      = /modules/apps/cplex/12.8/cplex
#CONCERTDIR    = /modules/apps/cplex/12.8/concert

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
C_EX =  raiz rf fo rfoks rfopr

all: $(C_EX) main
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o fo.o rfoks.o rfopr.o main.o -o HeurX $(CPLEX_FLAGS)
	rm -r *.o
	clear
	qsub script1
	qsub script2
	qsub script3
	qsub script3b
	qsub script4
	qsub script4b
	qsub script5
	qsub script6
	qsub script6b
	qsub script7
	qsub script8
	qsub script9
	qsub script10
	qsub script10b
	

# ------------------------------------------------------------
# The examples
#

raiz: raiz.h raiz.cpp
	$(CCC) -c $(CFLAGS) raiz.cpp -o raiz.o $(CPLEX_FLAGS)

rf: rf.h rf.cpp
	$(CCC) -c $(CFLAGS) rf.cpp -o rf.o $(CPLEX_FLAGS)

fo: fo.h fo.cpp
	$(CCC) -c $(CFLAGS) fo.cpp -o fo.o $(CPLEX_FLAGS)

rfoks: rfoks.h rfoks.cpp
	$(CCC) -c $(CFLAGS) rfoks.cpp -o rfoks.o $(CPLEX_FLAGS)

rfopr: rfopr.h rfopr.cpp
	$(CCC) -c $(CFLAGS) rfopr.cpp -o rfopr.o $(CPLEX_FLAGS)

main: main.cpp
	$(CCC) -c $(CFLAGS) main.cpp -o main.o $(CPLEX_FLAGS)
clean :
	rm -r *.o

# Local Variables:
# mode: makefile
# End:

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /opt/ibm/ILOG/CPLEX_Studio1271/cplex
CONCERTDIR    = /opt/ibm/ILOG/CPLEX_Studio1271/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0
CC  = gcc -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCOPT = -m64 -g -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC -fno-strict-aliasing
# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CLNFLAGS  = -lilocplex -lconcert -lcplex -lm -lpthread -fexceptions

.PHONY: all clean

#all:
#	make

#execute: all
#	make execute_c

CPLEXINCDIR   = $(CPLEXDIR)/include
EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
CFLAGS  = $(CCOPT)  -I$(CPLEXINCDIR)

#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX = raiz

all: $(C_EX) main_v2 
	$(CCC) raiz.o main_v2.o -o exec2 $(CLNFLAGS)
	./exec2

	 
# ------------------------------------------------------------
# The examples
#

raiz:raiz.h raiz.cpp
	$(CCC) -c $(CFLAGS) raiz.cpp -o raiz.o $(CLNFLAGS)

main_v2:main_v2.cpp
	$(CCC) -c $(CFLAGS) main_v2.cpp -o main_v2.o $(CLNFLAGS)

clean : 
	rm -r main.o

# Local Variables:
# mode: makefile
# End:

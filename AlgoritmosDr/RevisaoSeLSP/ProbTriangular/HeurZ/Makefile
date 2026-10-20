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
CPLEX_FLAGS = -lilocplex -lconcert -lcplex -lm -lpthread -fexceptions -lrt -ldl

.PHONY: all clean

CPLEXINCDIR   = $(CPLEXDIR)/include
CFLAGS  = $(COPT)  
#-I$(CPLEXINCDIR)

#------------------------------------------------------------
#  make all      : to compile the examples.
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX =  raiz rf fo rf2 fo2 rfoks rfopr

all: $(C_EX) main
	$(CCC) $(CFLAGS) raiz.o rf.o fo.o rf2.o fo2.o rfoks.o rfopr.o main.o -o allexec $(CPLEX_FLAGS)
	rm -r *.o
	clear
	#./allexec Data2-15-5-0.8-0.5-50-80-20-1-mod
	./allexec Data3-15-10-0.8-0.5-100-60-20-5-mod

# ------------------------------------------------------------
# The examples
#

raiz: raiz.h raiz.cpp
	$(CCC) -c $(CFLAGS) raiz.cpp -o raiz.o $(CPLEX_FLAGS)

rf: rf.h rf.cpp
	$(CCC) -c $(CFLAGS) rf.cpp -o rf.o $(CPLEX_FLAGS)

rf2: rfII.h rfII.cpp
	$(CCC) -c $(CFLAGS) rfII.cpp -o rf2.o $(CPLEX_FLAGS)

fo: fo.h fo.cpp
	$(CCC) -c $(CFLAGS) fo.cpp -o fo.o $(CPLEX_FLAGS)

fo2: foII.h foII.cpp
	$(CCC) -c $(CFLAGS) foII.cpp -o fo2.o $(CPLEX_FLAGS)

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

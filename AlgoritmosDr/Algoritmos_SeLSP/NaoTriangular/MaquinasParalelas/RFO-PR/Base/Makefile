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

RAIZ = /home/desymcb/Desiree/PIDLS/Raiz
BASE = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/Base
MAIN1 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-1
MAIN2 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-2
MAIN3 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-3
MAIN4 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-4
MAIN5 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-5
MAIN6 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-6
MAIN7 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-7
MAIN8 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-8
MAIN9 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-9
MAIN10 = /home/desymcb/Desiree/PIDLS/RFO-KS/KS-nb/RFO-10
#------------------------------------------------------------
#  make all      : to compile the examples.
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------
C_EX =  $(RAIZ)/raiz $(BASE)/rf $(MAIN1)/main $(MAIN2)/main $(MAIN3)/main $(MAIN4)/main $(MAIN5)/main $(MAIN6)/main $(MAIN7)/main $(MAIN8)/main $(MAIN9)/main $(MAIN10)/main

all: $(C_EX)
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN1)/main.o -o $(MAIN1)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN1)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN2)/main.o -o $(MAIN2)/rfo $(CPLEX_FLAGS)	
	clear
	qsub $(MAIN2)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN3)/main.o -o $(MAIN3)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN3)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN4)/main.o -o $(MAIN4)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN4)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN5)/main.o -o $(MAIN5)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN5)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN6)/main.o -o $(MAIN6)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN6)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN7)/main.o -o $(MAIN7)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN7)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN8)/main.o -o $(MAIN8)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN8)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN9)/main.o -o $(MAIN9)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN9)/script
	$(CCC) $(CFLAGS) $(CCLNDIRS) raiz.o rf.o $(MAIN10)/main.o -o $(MAIN10)/rfo $(CPLEX_FLAGS)
	clear
	qsub $(MAIN10)/script

# ------------------------------------------------------------
# The examples
#

$(RAIZ)/raiz: $(RAIZ)/raiz.h $(RAIZ)/raiz.cpp
	$(CCC) -c $(CFLAGS) $(RAIZ)/raiz.cpp -o raiz.o $(CPLEX_FLAGS)

$(BASE)/rf: $(BASE)/rf.h $(BASE)/rf.cpp
	$(CCC) -c $(CFLAGS) $(BASE)/rf.cpp -o rf.o $(CPLEX_FLAGS)

$(MAIN1)/main: $(MAIN1)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN1)/main.cpp -o $(MAIN1)/main.o $(CPLEX_FLAGS)

$(MAIN2)/main: $(MAIN2)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN2)/main.cpp -o $(MAIN2)/main.o $(CPLEX_FLAGS)

$(MAIN3)/main: $(MAIN3)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN3)/main.cpp -o $(MAIN3)/main.o $(CPLEX_FLAGS)

$(MAIN4)/main: $(MAIN4)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN4)/main.cpp -o $(MAIN4)/main.o $(CPLEX_FLAGS)

$(MAIN5)/main: $(MAIN5)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN5)/main.cpp -o $(MAIN5)/main.o $(CPLEX_FLAGS)

$(MAIN6)/main: $(MAIN6)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN6)/main.cpp -o $(MAIN6)/main.o $(CPLEX_FLAGS)

$(MAIN7)/main: $(MAIN7)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN7)/main.cpp -o $(MAIN7)/main.o $(CPLEX_FLAGS)

$(MAIN8)/main: $(MAIN8)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN8)/main.cpp -o $(MAIN8)/main.o $(CPLEX_FLAGS)

$(MAIN9)/main: $(MAIN9)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN9)/main.cpp -o $(MAIN9)/main.o $(CPLEX_FLAGS)

$(MAIN10)/main: $(MAIN10)/main.cpp
	$(CCC) -c $(CFLAGS) $(MAIN10)/main.cpp -o $(MAIN10)/main.o $(CPLEX_FLAGS)
clean :
	rm -r *.o

# Local Variables:
# mode: makefile
# End:

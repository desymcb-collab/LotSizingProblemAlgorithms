#ifndef RAIZ_H_INCLUDED
#define RAIZ_H_INCLUDED

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <limits.h>

#define	INF    999999999.0
#define zero        0.00001
#define Max_T   16
#define Max_M   3
#define BigM    25
#define vetKS 15

#define LINE_BUF_LEN     255
#define PROG_ID_STR      "\nRFO algorithm!!\n"

#ifndef STR_ERR_UNKNOWN_LONG_OPT
# define STR_ERR_UNKNOWN_LONG_OPT   "%s: unrecognized option `--%s'\n"
#endif

#ifndef STR_ERR_LONG_OPT_AMBIGUOUS
# define STR_ERR_LONG_OPT_AMBIGUOUS "%s: option `--%s' is ambiguous\n"
#endif

#ifndef STR_ERR_MISSING_ARG_LONG
# define STR_ERR_MISSING_ARG_LONG   "%s: option `--%s' requires an argument\n"
#endif

#ifndef STR_ERR_UNEXPEC_ARG_LONG
# define STR_ERR_UNEXPEC_ARG_LONG   "%s: option `--%s' doesn't allow an argument\n"
#endif

#ifndef STR_ERR_UNKNOWN_SHORT_OPT
# define STR_ERR_UNKNOWN_SHORT_OPT  "%s: unrecognized option `-%c'\n"
#endif

#ifndef STR_ERR_MISSING_ARG_SHORT
# define STR_ERR_MISSING_ARG_SHORT  "%s: option `-%c' requires an argument\n"
#endif

#define STR_HELP_JUMP \
  "  -a, --jump          # jump (qtd de periodos movimentados no RF)\n"

#define STR_HELP_JUMPFO \
  "  -b, --jumpfo           # jumpfo (qtd de periodos movimentados no FO)\n"

#define STR_HELP_INT \
  "  -c, --binary            # binary: qtd de periodos com dominio inteiro\n"

#define STR_HELP_QUIET \
  "      --quiet           reduce output to a minimum, no extra files\n"

#define STR_HELP_HELP \
  "  -h, --help            display this help text and exit\n"

#define STR_HELP_TSPLIBFILE \
  "  -i, --tsplibfile     f inputfile (TSPLIB format necessary)\n"

  #define STR_HELP_SEED \
  "      --seed           # seed for the random number generator\n"

static const char *const STR_HELP[] = {
  STR_HELP_JUMP ,
  STR_HELP_JUMPFO ,
  STR_HELP_INT ,
  STR_HELP_QUIET ,
  STR_HELP_HELP,
  STR_HELP_TSPLIBFILE ,
  STR_HELP_SEED       ,
  NULL
};

struct options {
/* Set to 1 if option --jump (-a) has been specified.  */
unsigned int opt_jump : 1;

/* Set to 1 if option --jumpfo (-b) has been specified.  */
unsigned int opt_jumpfo : 1;

/* Set to 1 if option --binary (-c) has been specified.  */
unsigned int opt_binary : 1;

/* Set to 1 if option --quiet has been specified.  */
unsigned int opt_quiet : 1;

/* Set to 1 if option --help (-h) has been specified.  */
unsigned int opt_help : 1;

/* Set to 1 if option --seed has been specified.  */
unsigned int opt_seed : 1;

/* Argument to option --jump (-a).  */
const char *arg_jump;

/* Argument to option --jumpfo (-b).  */
const char *arg_jumpfo;

/* Argument to option --binary (-c).  */
const char *arg_binary;

/* Argument to option --tsplibfile (-i).  */
const char *arg_tsplibfile;

/* Set to 1 if option --tsplibfile (-i) has been specified.  */
unsigned int opt_tsplibfile : 1;

/* Argument to option --seed.  */
const char *arg_seed;

};

static int parse_options (struct options *const options,const char *const program_name, const int argc, char **const argv);
static void check_out_of_range (double value, double minval, double maxval, const char *optionName);

typedef struct guarda{
    double sol;
    std::vector< std::vector< std::vector < double > > > xstore;
	std::vector< std::vector< std::vector < int > > > ystore;
	std::vector< std::vector< double > > zstore;
	//std::vector< std::vector< std::vector<double> > > istore;
} guarda;

using namespace std;

class Instance{
    public:
        Instance(){}
        ~Instance(){}

    /************************************************************************
                                   Parameters
    ************************************************************************/
    ///Paralell machine
    std::vector< double > Capacity;                      //plant capacity MxT
    std::vector< double > HoldingCost;                                  //holding cost N
    std::vector< std::vector< double > > SetupTime;      //setup time NxNxM
    std::vector< std::vector< double > > SetupCost;      //setup cost NxNxM
    std::vector< double > ProdTime;                      //production time NxM
    std::vector< double > ProdCost;                      //production time NxM
    std::vector< std::vector< double > > Demand;                        //production demand NxT
    std::vector< std::vector< double > > MaxProd;                       //production cost NxT
    std::vector< int > MachineProd;                   //binary variable - 1 indicates that machine k produced item i, 0 otherwise MxN
    std::vector< int> isetup;    //estado inicial da maquina para cada produto
    std::vector< int> istock;    //estoque inicial de cada produto
    std::vector< int > minl;      //tamanho minimo do lote
    std::vector< double> maxl;      //tamanho maximo do lote
    // int   			p,p1;   //parameter from feasibility phase
    int             T, //number of periods
                    M, //number of plants
                    N; //number of itens

    double          solfea; //feasibility process solution
    int             jump, jump_fo,binary;

    char            filename[80]; //name of the output solutions
    char            pathname[300]; //name of the path in which the filename will be storaged
    char            fileinfo[300]; //name of the file with the instances information

    FILE *info;
    FILE *info2;

    /************************************************************************
                                    Functions
    ************************************************************************/

    void inicializa_variaveis(void);
    void desaloca_variaveis(void);
    void le_dados_SM(char * input_data);
    void le_dados_M(char * input_data);
    void show_instance(void);

    void set_default_parameters(void);
    int parse_commandline (int argc, char *argv []);
    void init_program(long int argc, char *argv[] );


};

#endif // RAIZ_H_INCLUDED

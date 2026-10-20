#ifndef RAIZ_H_INCLUDED
#define RAIZ_H_INCLUDED

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#define INF     999999999.0
#define zero    0.00001
#define Max_N   26
#define Max_T   16
#define BigM    25
#define vetKS   7

typedef struct guarda{
    double sol;
    std::vector< std::vector < std::vector<double> > > xstore;
    std::vector< std::vector< std::vector <double> > > ystore;
    std::vector< std::vector<double> > zstore;
    //std::vector< std::vector< std::vector<double> > > Istore;
} guarda;

/************************************************************************
                    Estutura de dados kernel search
************************************************************************/

struct dado_ks {
    int i1; //item atual
    int t1; //indica qual caminho foi escolhido (periodo)
    int j1; //indica qual caminho foi escolhido (item)
    double cost; //fornece o valor deste caminho
    int kernel;
};

struct LRC_ks {
    int nelementos; //indica o número de elementos da lrc
    dado_ks campo[Max_N * Max_N * Max_T]; //dados da estrutura
};

using namespace std;

class Instance{
    public:
        Instance(){}
        ~Instance(){}

    /************************************************************************
                                   Variables
    ************************************************************************/

    std::vector< std::vector< std::vector <double> > > X;
    std::vector< std::vector< std::vector <double> > > Y;
    std::vector< std::vector< std::vector <double> > > reducedCostY;
    std::vector< std::vector<double> > Z;
    std::vector< std::vector<double> > G;
    std::vector< std::vector< std::vector <double> > > F;

    /************************************************************************
                                   Parameters
    ************************************************************************/
    ///Paralell machine
    std::vector< double > Capacity;                                     //plant capacity T
    std::vector< double > HoldingCost;                                  //holding cost N
    std::vector< std::vector< double > > SetupTime;                     //setup time NxN
    std::vector< std::vector< double > > SetupCost;                     //setup cost NxN
    std::vector< double > ProdTime;                                     //production time N
    std::vector< double > ProdCost;                                     //production time N
    std::vector< std::vector< double > > Demand;                        //production demand NxT
    std::vector< std::vector< double > > MaxProd;                       //production cost NxT
    std::vector< int > MachineProd;                      //binary variable - 1 indicates that machine k produced item i, 0 otherwise N
    std::vector< double> isetup;                                        //estado inicial da maquina para cada produto
    std::vector< double> istock;                                        //estoque inicial de cada produto
    std::vector< int > minl;                             //tamanho minimo do lote
    std::vector< double> maxl;                                          //tamanho maximo do lote
    
    int             T,  //number of periods
                    N,  //number of itens
                    nontri[3];  //nontriangular items

    double          solutionRF,
                    solutionRFO,
                    solutionKS,
                    solutionRFKS,
                    solutionKSPR,
                    solutionRFPR,
                    LowerBound;

    char            filename[80]; //name of the output solutions
    char            pathname[80]; //name of the path in which the filename will be storaged
    char            fileinfo[200]; //name of the file with the instances information

    FILE *info;
    FILE *info2;
    guarda vet[vetKS];          //keep the linear relaxation in vet[0] and the best solution in vet[1]

    /************************************************************************
                                    Functions
    ************************************************************************/

    void inicializa_variaveis(void);
    void desaloca_variaveis(void);
    void le_dados_SM(const char * input_data);
    void show_instance(void);
    void show_solution(void);
};

#endif // RAIZ_H_INCLUDED
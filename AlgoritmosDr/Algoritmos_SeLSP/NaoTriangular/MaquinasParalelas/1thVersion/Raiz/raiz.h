#ifndef RAIZ_H_INCLUDED
#define RAIZ_H_INCLUDED

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#define	INF    999999999.0
#define zero        0.00001
#define Max_T   16
#define Max_M   3
#define BigM    25
#define vetKS 15

typedef struct guarda{
    double sol;
    std::vector< std::vector< std::vector < std::vector<double> > > > xstore;
	std::vector< std::vector< std::vector < std::vector<int> > > > ystore;
	std::vector< std::vector< std::vector<double> > > zstore;
	//std::vector< std::vector< std::vector<double> > > Istore;
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
    std::vector< std::vector< double > > Capacity;                      //plant capacity MxT
    std::vector< double > HoldingCost;                                  //holding cost N
    std::vector< std::vector< std::vector< double > > > SetupTime;      //setup time NxNxM
    std::vector< std::vector< std::vector< double > > > SetupCost;      //setup cost NxNxM
    std::vector< std::vector< double > > ProdTime;                      //production time NxM
    std::vector< std::vector< double > > ProdCost;                      //production time NxM
    std::vector< std::vector< double > > Demand;                        //production demand NxT
    std::vector< std::vector< double > > MaxProd;                       //production cost NxT
    std::vector< std::vector< int > > MachineProd;                   //binary variable - 1 indicates that machine k produced item i, 0 otherwise MxN
    std::vector< double> isetup;    //estado inicial da maquina para cada produto
    std::vector< double> istock;    //estoque inicial de cada produto
    std::vector< std::vector< int > > minl;      //tamanho minimo do lote
    std::vector< double> maxl;      //tamanho maximo do lote
    // int   			p,p1;   //parameter from feasibility phase
    int             T, //number of periods
                    M, //number of plants
                    N; //number of itens

    double          solfea; //feasibility process solution

    char            filename[80]; //name of the output solutions
    char            pathname[80]; //name of the path in which the filename will be storaged
    char            fileinfo[200]; //name of the file with the instances information

    FILE *info;
    FILE *info2;

    /************************************************************************
                                    Functions
    ************************************************************************/

    void inicializa_variaveis(void);
    void desaloca_variaveis(void);
    void le_dados_SM(char * input_data);
    void show_instance(void);
};

#endif // RAIZ_H_INCLUDED

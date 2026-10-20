#include "raiz.h"

/************Inicializa variaveis do problema************/
void Instance::inicializa_variaveis(void)
{
    Capacity.resize(T+1,0.0);

    HoldingCost.resize(N+1,0);
    SetupTime.resize(N+1);
    SetupCost.resize(N+1);
    ProdCost.resize(N+1,0.0);
    ProdTime.resize(N+1,0.0);
    Demand.resize(N+1);
    MaxProd.resize(N+1);
    MachineProd.resize(N+1,0.0);
    minl.resize(N+1,0.0);
    isetup.resize(N+1,0);
    istock.resize(N+1,0);
    for(int i=0;i<=N;i++)
    {
        Demand[i].resize(T+1,0.0);
        MaxProd[i].resize(T+2,0.0);
        SetupTime[i].resize(N+1,0.0);
        SetupCost[i].resize(N+1,0.0);
    }

}

/************Inicializa variaveis do problema************/
void Instance::desaloca_variaveis(void)
{
    for(int i=0;i<=N;i++)
    {
        Demand[i].clear();
        MaxProd[i].clear();
        SetupTime[i].clear();
        SetupCost[i].clear();
    }

    Capacity.clear();
    HoldingCost.clear();
    SetupTime.clear();
    SetupCost.clear();
    ProdCost.clear();
    ProdTime.clear();
    Demand.clear();
    MaxProd.clear();
    MachineProd.clear();
    minl.clear();
    isetup.clear();
    istock.clear();
}

/********************Read instances**********************/
void Instance::le_dados_SM(char * input_data)
{
    //printf("entrou\n");
    FILE   *arq_ent;

    double  faux1, faux2, faux3, faux4;
    double  aux1;
    int Iaux;

    arq_ent = fopen(input_data,"r");

    //printf("leu\n");
    int r3;
    char aux[20];
    /* read general data of the problem */
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&N);
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&T);

    inicializa_variaveis();
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        HoldingCost[i] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);

    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        ProdTime[i] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);

    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            SetupTime[i][j] = aux1;
        }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
        for(int j=1; j<=N; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            SetupCost[i][j] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);

    for(int t=1; t<=T; t++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        Capacity[t] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            Demand[i][t] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=T; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            MaxProd[i][j] = aux1;
        }
        MaxProd[i][T+1] = 0;
    }

    r3 = fscanf(arq_ent,"%s",aux);

    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%d",&Iaux);
        isetup[i] = Iaux;
        MachineProd[i] = 1;
        //printf("%d ",MachineProd[i][k]);
    }
        //printf("\n");
    r3 = fscanf(arq_ent,"%s",aux);

    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%d",&Iaux);
        istock[i] = Iaux;
        //printf("%d ",MachineProd[i][k]);
    }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%d",&Iaux);
        minl[i] = Iaux;
    }

    fclose(arq_ent);
}

/*void Instance::le_dados_M(char * input_data)
{
    FILE   *arq_ent;

    double  faux1, faux2, faux3, faux4;
    double  aux1;
    int Iaux;
    //printf("Entrou!!\n");

    arq_ent = fopen(input_data,"r");

    int r3;
    char aux[20];
     read general data of the problem
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&N);
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&T);
    r3 = fscanf(arq_ent,"%s",aux);
    r3 = fscanf(arq_ent,"%d",&M);

    inicializa_variaveis();
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
    {
        r3 = fscanf(arq_ent,"%lf",&aux1);
        HoldingCost[i] = aux1;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            ProdTime[i][k] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
            {
                r3 = fscanf(arq_ent,"%lf",&aux1);
                SetupTime[i][j][k] = aux1;
            }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int j=1; j<=N; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            SetupTime[0][j][k] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            Capacity[k][t] = aux1;
        }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int i=1; i<=N; i++)
        for(int t=1; t<=T; t++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            Demand[i][t] = aux1;
        }

    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
        {
            MaxProd[i][t] = 0;
            for(int u=1; u<=T; u++)
                MaxProd[i][t] += Demand[i][u];
        }
        MaxProd[i][T+1] = 0;
    }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
            for(int j=1; j<=N; j++)
            {
                r3 = fscanf(arq_ent,"%lf",&aux1);
                SetupCost[i][j][k] = aux1;
            }
    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int j=1; j<=N; j++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            SetupCost[0][j][k] = aux1;
        }
    //r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
    {
        for(int i=1; i<=N; i++)
        {
            //r3 = fscanf(arq_ent,"%d",&Iaux);
            MachineProd[i][k] = 1;
            //printf("%d ",MachineProd[i][k]);
        }
        //printf("\n");
    }

    r3 = fscanf(arq_ent,"%s",aux);
    for(int k=1; k<=M; k++)
        for(int i=1; i<=N; i++)
        {
            r3 = fscanf(arq_ent,"%lf",&aux1);
            minl[i][k] = aux1;
        }

    fclose(arq_ent);
}*/

/**************Mostra dados das instancias***************/
void Instance::show_instance(void)
{
    printf("%d\n",N);
    printf("%d\n",T);


    for(int i=1; i<=N; i++)
        printf("%.1lf ",HoldingCost[i]);

    printf("\n");

        for(int i=1; i<=N; i++)
            printf("%.1lf ",ProdTime[i]);
        printf("\n");


    printf("\n");


        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
                printf("%.1lf ",SetupTime[i][j]);
            printf("\n");
        }

        for(int t=1; t<=T; t++)
            printf("%.1lf ",Capacity[t]);
        printf("\n");


    printf("\n");
    for(int i=1; i<=N; i++)
    {
        for(int t=1; t<=T; t++)
            printf("%.1lf ",Demand[i][t] );
        printf("\n");
    }

    for(int i=1; i<=N; i++)
    {
        for(int j=1; j<=T; j++)
            printf("%.1lf ",MaxProd[i][j]);
        printf("\n");
    }

        for(int i=1; i<=N; i++)
        {
            for(int j=1; j<=N; j++)
                printf("%.1lf ",SetupCost[i][j]);
            printf("\n");
        }

        for(int i=1; i<=N; i++)
        {
            printf("%d ",MachineProd[i]);
        }
        printf("\n");

}


void Instance:: set_default_parameters(void)
/*
      FUNCTION: set default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    jump = 1;
    jump_fo = 1;
    binary = 3;
}


/* Parse command line options.  Return index of first non-option argument,
   or -1 if an error is encountered.  */
static int parse_options (struct options *const options,const char *const program_name,const int argc, char **const argv)
{
  static const char *const optstr__jump = "jump";
  static const char *const optstr__jumpfo = "jumpfo";
  static const char *const optstr__binary = "binary";
  static const char *const optstr__quiet = "quiet";
  static const char *const optstr__help = "help";
  static const char *const optstr__tsplibfile = "tsplibfile";
  static const char *const optstr__seed = "seed";
  int i = 0;
  options->opt_jump = 0;
  options->opt_jumpfo = 0;
  options->opt_binary = 0;
  options->opt_quiet = 0;
  options->opt_help = 0;
  options->arg_jump = 0;
  options->arg_jumpfo = 0;
  options->arg_binary = 0;
  options->arg_tsplibfile = 0;
  options->opt_seed = 0;
  options->arg_seed = 0;
  while (++i < argc)
  {
    const char *option = argv [i];
    if (*option != '-')
      return i;
    else if (*++option == '\0')
      return i;
    else if (*option == '-')
    {
      const char *argument;
      size_t option_len;
      ++option;
      if ((argument = strchr (option, '=')) == option)
        goto error_unknown_long_opt;
      else if (argument == 0)
        option_len = strlen (option);
      else
        option_len = argument++ - option;
      switch (*option)
      {
       case '\0':
        return i + 1;
       case 'a':
        if (strncmp (option + 1, optstr__jump + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
            options->arg_jump = argument;
          else if (++i < argc)
            options->arg_jump = argv [i];
          else
          {
            option = optstr__jump;
            goto error_missing_arg_long;
          }
          options->opt_jump = 1;
          break;
        }
        goto error_unknown_long_opt;
       case 'b':
        if (strncmp (option + 1, optstr__jumpfo + 1, option_len - 1) == 0)
        {
          if (option_len <= 1)
            goto error_long_opt_ambiguous;
          if (argument != 0)
            options->arg_jumpfo = argument;
          else if (++i < argc)
            options->arg_jumpfo = argv [i];
          else
          {
            option = optstr__jumpfo;
            goto error_missing_arg_long;
          }
          options->opt_jumpfo = 1;
          break;
        }
        goto error_unknown_long_opt;
       case 'c':
        if (strncmp (option + 1, optstr__binary + 1, option_len - 1) == 0)
        {
          if (argument != 0)
            options->arg_binary = argument;
          else if (++i < argc)
            options->arg_binary = argv [i];
          else
          {
            option = optstr__binary;
            goto error_missing_arg_long;
          }
          options->opt_binary = 1;
          break;
        }
        goto error_unknown_long_opt;
       case 'h':
        if (strncmp (option + 1, optstr__help + 1, option_len - 1) == 0)
        {
          if (argument != 0)
          {
            option = optstr__help;
            goto error_unexpec_arg_long;
          }
          options->opt_help = 1;
          return i + 1;
        }
        goto error_unknown_long_opt;
       case 'q':
        if (strncmp (option + 1, optstr__quiet + 1, option_len - 1) == 0)
        {
            if (argument != 0)
            {
                option = optstr__quiet;
                goto error_unexpec_arg_long;
            }
            options->opt_quiet = 1;
            break;
        }
        goto error_unknown_long_opt;
        case 's':
          if (strncmp (option + 1, optstr__seed + 1, option_len - 1) == 0)
          {
              if (option_len <= 1)
                  goto error_long_opt_ambiguous;
              if (argument != 0)
                  options->arg_seed = argument;
              else if (++i < argc)
                  options->arg_seed = argv [i];
              else
              {
                  option = optstr__seed;
                  goto error_missing_arg_long;
              }
              options->opt_seed = 1;
              break;
          }
          goto error_unknown_long_opt;
       default:
       error_unknown_long_opt:
        fprintf (stderr, STR_ERR_UNKNOWN_LONG_OPT, program_name, option);
        return -1;
       error_long_opt_ambiguous:
        fprintf (stderr, STR_ERR_LONG_OPT_AMBIGUOUS, program_name, option);
        return -1;
       error_missing_arg_long:
        fprintf (stderr, STR_ERR_MISSING_ARG_LONG, program_name, option);
        return -1;
       error_unexpec_arg_long:
        fprintf (stderr, STR_ERR_UNEXPEC_ARG_LONG, program_name, option);
        return -1;
      }
    }
    else
      do
      {
        switch (*option)
        {
         case 'a':
          if (option [1] != '\0')
            options->arg_jump = option + 1;
          else if (++i < argc)
            options->arg_jump = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_jump = 1;
          break;
         case 'b':
          if (option [1] != '\0')
            options->arg_jumpfo = option + 1;
          else if (++i < argc)
            options->arg_jumpfo = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_jumpfo = 1;
          break;
         case 'c':
          if (option [1] != '\0')
            options->arg_binary = option + 1;
          else if (++i < argc)
            options->arg_binary = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_binary = 1;
          break;
          case 'i':
          if (option [1] != '\0')
            options->arg_tsplibfile = option + 1;
          else if (++i < argc)
            options->arg_tsplibfile = argv [i];
          else
            goto error_missing_arg_short;
          option = "\0";
          options->opt_tsplibfile = 1;
          break;
         default:
         fprintf (stderr, STR_ERR_UNKNOWN_SHORT_OPT, program_name, *option);
          return -1;
         error_missing_arg_short:
          fprintf (stderr, STR_ERR_MISSING_ARG_SHORT, program_name, *option);
          return -1;
        }
      } while (*++option != '\0');
  }
  return i;
}

static void check_out_of_range (double value, double minval, double maxval,
                    const char *optionName)
/*
      FUNCTION: check whether parameter values are within allowed range
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    if (value < minval || value > maxval) {
        fprintf (stderr,"Error: Option `%s' out of range [%g, %g]\n",
                 optionName, minval, maxval);
        exit(1);
    }
}

int Instance :: parse_commandline (int argc, char *argv [])
{
    int i;
    const char *progname;
    struct options options;

    progname = argv [0] != NULL && *(argv [0]) != '\0'
	? argv [0]
	: "rfo";

    i = parse_options (&options,progname, argc, argv);

    if (i < 2)
    {
	fprintf (stderr, "No options are specified\n");
	fprintf (stderr, "Try `%s --help' for more information.\n",
		 progname);
	exit(1);
    }

    if (options.opt_help)
    {
        int k;
	printf ("Usage: %s [OPTION]... [ARGUMENT]...\n"
		"Options:\n", progname);
        for (k = 0; STR_HELP[k]; k++)
            printf ("%s", STR_HELP[k]);
	exit(0);
    }

    puts ("\t OPTIONS:");


    if ( options.opt_jump ) {
	jump = int(atof(options.arg_jump));
	fputs ("  -a  --jump ", stdout);
	if (options.arg_jump != NULL)
	    printf ("with argument \"%d\"\n", jump);
	check_out_of_range( jump, 1., 5., "jump");
    } else {
        fprintf(stdout,"\tNote: jump is set to default %d\n", jump);
    }

    if ( options.opt_jumpfo ) {
	jump_fo = int(atof(options.arg_jumpfo));
	fputs ("  -b  --jumpfo ", stdout);
	if (options.arg_jumpfo != NULL)
	    printf ("with argument \"%d\"\n", jump_fo);
	check_out_of_range( jump_fo, 1., 5., "jumpfo");
    } else {
        fprintf(stdout,"\tNote: jumpfo is set to default %d\n", jump_fo);
    }


    if ( options.opt_binary ) {
	binary = int(atof(options.arg_binary));
	fputs ("  -c  --binary ", stdout);
	if (options.arg_binary != NULL)
	    printf ("with argument \"%d\"\n", binary);
        check_out_of_range( binary, 1., 5., "binary");
    } else {
        fprintf(stdout,"\tNote: binary is set to default %d\n", binary);
    }

//    if ( options.opt_seed ) {
//        seed = atol(options.arg_seed);
//        fputs ("      --seed ", stdout);
//        if (options.arg_seed != NULL)
//            printf ("with argument \"%ld\"\n", seed);
//    } else {
//        fprintf(stdout,"\tNote: a seed was generated as %ld\n", seed);
//    }

    if ( options.opt_tsplibfile )
    {
        if (strlen(options.arg_tsplibfile) >= LINE_BUF_LEN) {
            fprintf (stderr, "error: too long input filename '%s', maximum length is %d", options.arg_tsplibfile, LINE_BUF_LEN);
            exit (1);
        }
	strcpy (fileinfo, options.arg_tsplibfile);
 	fputs ("  -i  --tsplibfile ", stdout);
	if (options.arg_tsplibfile != NULL)
	    printf ("with argument \"%s\"\n", fileinfo );
    }

    puts ("Non-option arguments:");

    while (i < argc) {
	fprintf (stderr,"  \"%s\"\n", argv [i++]);
	fprintf (stderr,"\nThere were non-option arguments\n");
	fprintf (stderr,"I suspect there is something wrong, maybe wrong option name; exit\n");
	exit(1);
    }

    return 0;
}




void Instance :: init_program(long int argc, char *argv[] )
/*
      FUNCTION:       initialize the program,
      INPUT:          program arguments, needed for parsing commandline
      OUTPUT:         none
      COMMENTS:
*/
{
  char temp_buffer[LINE_BUF_LEN];

  printf(PROG_ID_STR);
  set_default_parameters();
  setbuf(stdout,NULL);
  parse_commandline(argc, argv);
  printf("\n.. done\n");
}

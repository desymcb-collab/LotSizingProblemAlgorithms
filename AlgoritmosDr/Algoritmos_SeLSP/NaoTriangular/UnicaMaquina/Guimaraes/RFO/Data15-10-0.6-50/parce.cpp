/* This file has been generated with opag 0.6.4.  */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include "ls.h"


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

static const char *const STR_HELP[] = {
  STR_HELP_JUMP ,
  STR_HELP_JUMPFO ,
  STR_HELP_INT ,
  STR_HELP_QUIET ,
  STR_HELP_HELP,
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

/* Argument to option --jump (-a).  */
const char *arg_jump;

/* Argument to option --jumpfo (-b).  */
const char *arg_jumpfo;

/* Argument to option --binary (-c).  */
const char *arg_binary;

};

/* Parse command line options.  Return index of first non-option argument,
   or -1 if an error is encountered.  */
static int
parse_options (struct options *const options, const char *const program_name,
               const int argc, char **const argv)
{
  static const char *const optstr__jump = "jump";
  static const char *const optstr__jumpfo = "jumpfo";
  static const char *const optstr__binary = "binary";
  static const char *const optstr__quiet = "quiet";
  static const char *const optstr__help = "help";
  int i = 0;
  options->opt_jump = 0;
  options->opt_jumpfo = 0;
  options->opt_binary = 0;
  options->opt_quiet = 0;
  options->opt_help = 0;
  options->arg_jump = 0;
  options->arg_jumpfo = 0;
  options->arg_binary = 0;
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

static void
check_out_of_range (double value, double minval, double maxval,
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

int parse_commandline (RF *rf ,int argc, char *argv [])
{
    int i;
    const char *progname;
    struct options options;

    progname = argv [0] != NULL && *(argv [0]) != '\0'
	? argv [0]
	: "rfo";

    i = parse_options (&options, progname, argc, argv);

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

    quiet_flag = options.opt_quiet;

    if (ls_flag) {
        set_default_ls_parameters();
    }


    if ( options.opt_jump ) {
	rf->jump = atof(options.arg_jump);
	fputs ("  -a  --jump ", stdout);
	if (options.arg_jump != NULL)
	    printf ("with argument \"%f\"\n", rf->jump);
	check_out_of_range( rf->jump, 0., 100., "jump");
    } else {
        fprintf(stdout,"\tNote: jump is set to default %g\n", rf->jump);
    }

    if ( options.opt_jumpfo ) {
	rf->jump_fo = atof(options.arg_jumpfo);
	fputs ("  -b  --jumpfo ", stdout);
	if (options.arg_jumpfo != NULL)
	    printf ("with argument \"%f\"\n", rf->jumpfo);
	check_out_of_range( rf->jumpfo, 0., 100., "jumpfo");
    } else {
        fprintf(stdout,"\tNote: jumpfo is set to default %g\n", rf->jumpfo);
    }


    if ( options.opt_binary ) {
	rf->binary = atof(options.arg_binary);
	fputs ("  -c  --binary ", stdout);
	if (options.arg_binary != NULL)
	    printf ("with argument \"%f\"\n", rf->binary);
        check_out_of_range( rf->binary, 0.000001, 1., "binary");
    } else {
        fprintf(stdout,"\tNote: binary is set to default %g\n", rf->binary);
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

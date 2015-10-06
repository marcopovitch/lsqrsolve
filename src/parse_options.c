#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <popt.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "parse_options.h"

void lsqrsolve_info()
{
    char *rayver = librayversion();
    char *meshver = libmeshversion();
    char *sparsever = libsparseversion();

    fprintf(stderr, "%s: version %s\n", PACKAGE, VERSION);
    fprintf(stderr, "using");
    fprintf(stderr, "\t%s\n", rayver);
    fprintf(stderr, "\t%s\n", meshver);
    fprintf(stderr, "\t%s\n", sparsever);

    free(rayver);
    free(meshver);
    free(sparsever);
}

/*
 * \brief parse command line using popt : checks validity of options
 * specified and returns parameters values.
 *
 * PRECOND : the default values for optional arguments must be set.
 */
void parse_command_line(int argc, char **argv,
                        char **mesh_filename,
                        char **vmodel,
                        char **importfilename,
                        char **log_filename,
                        char **output_filename,
                        char **output_type,
                        int *max_iter,
                        double *damping,
                        double *grad_damping,
                        int *use_ach, int *check_sparse)
{

    /* these intermediate var are only */
    /* there to build on IRIX with popt */
    static char *meshfile_opt = NULL;
    static char *vmodel_opt = NULL;
    static char *outputfile_opt = NULL;
    static char *importfilename_opt = NULL;
    static char *outputtype_opt = "gms";        /* default value */
    static char *logfile_opt = NULL;

    static int max_iter_opt = -1;
    static float damping_opt = 0.0;
    static float grad_damping_opt = 0.0;

    int rc;
    poptContext cx;

#include "options.h"

    /* args parsing using popt */
    cx = poptGetContext(PACKAGE, argc, (const char **) argv, options, 0);
    if (argc == 1) {
        poptPrintUsage(cx, stdout, 0);
        exit(1);
    }

    do {
        rc = poptGetNextOpt(cx);
        switch (rc) {
        case OPT_HELP:
            poptPrintHelp(cx, stdout, 0);
            exit(0);
        case OPT_VERSION:
            lsqrsolve_info();
            exit(0);
        case OPT_VELOCITY_MODEL:
            if (access(vmodel_opt, R_OK) < 0) {
                fprintf(stderr, "Can't access to velocity file '%s'\n",
                        vmodel_opt);
                exit(1);
            }
            break;
        case OPT_OUTPUT_TYPE:
            if (!strchr(outputtype_opt, 'm') &&
                !strchr(outputtype_opt, 'g') &&
                !strchr(outputtype_opt, 's')) {
                fprintf(stderr,
                        "output_type must be 'm'atlab, 'g'mt or 's'co.\n");
                exit(1);
            }
            break;
        case OPT_LOG_FILE:
            if (!logfile_opt) {
                fprintf(stderr,
                        "a log filename must be provided avec -l.\n");
                exit(1);
            }
            if (logfile_opt[0] == '-') {
                fprintf(stderr,
                        "Logfilenames may not begin with a - character.\n");
                exit(1);
            }
            break;
        case OPT_DAMPING:
            /*fprintf(stderr, "damping %f\n", damping_opt); */
            break;
        case OPT_GRAD_DAMPING:
            /*fprintf(stderr, "gradient damping %f\n", grad_damping_opt); */
            break;
        case OPT_CHECK_SPARSE:
            *check_sparse = 1;
            break;
        case OPT_USE_ACH:
            *use_ach = 1;
            break;
        case POPT_ERROR_BADOPT:        /* error */
            fprintf(stderr, "%s: %s\n", poptBadOption(cx, 0),
                    poptStrerror(rc));
            poptPrintUsage(cx, stdout, 0);
            exit(1);
        }
    } while (rc > 0);

    if (!meshfile_opt || !importfilename_opt || !outputfile_opt) {
        fprintf(stderr, "-m, -i and -o  must be provided\n");
        exit(1);
    }

    *mesh_filename = meshfile_opt;
    *vmodel = vmodel_opt;
    *importfilename = importfilename_opt;
    *output_filename = outputfile_opt;
    *log_filename = logfile_opt;
    *output_type = outputtype_opt;

    *damping = damping_opt;
    *grad_damping = grad_damping_opt;
    *max_iter = max_iter_opt;
}

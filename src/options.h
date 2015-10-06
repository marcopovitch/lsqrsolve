enum options {
    OPT_ZERO,                   /* must be 0 */
    OPT_HELP,
    OPT_MESHFILE,
    OPT_VELOCITY_MODEL,
    OPT_OUTPUT_FILE,
    OPT_OUTPUT_TYPE,
    OPT_USE_ACH,
    OPT_CHECK_SPARSE,
    OPT_IMPORTFILE,
    OPT_LOG_FILE,
    OPT_MAX_ITER,
    OPT_DAMPING,
    OPT_GRAD_DAMPING,
    OPT_VERSION
};

const struct poptOption options[] = {
    /* long optin,short option, type,var,default val,explanation,expl 2 */
    {
     "help", 'h', POPT_ARG_NONE,
     NULL, OPT_HELP,
     "help", "display this help message"},
    {
     "meshfile", 'm', POPT_ARG_STRING,
     &meshfile_opt, OPT_MESHFILE,
     "mesh config file (xml)", "FILE"},
    {
     "vmodel", '\0', POPT_ARG_STRING,
     &vmodel_opt, OPT_VELOCITY_MODEL,
     "velocity model", "FILE"},
    {
     "format", 'f', POPT_ARG_STRING,
     &outputtype_opt, OPT_OUTPUT_TYPE,
     "output format (matlab, gmt, sco)", "m|g|s"},
    {
     "ach", 'a', POPT_ARG_NONE,
     NULL, OPT_USE_ACH,
     "seismic ACH inversion", "use ACH"},
    {
     "check-sparse", 'c', POPT_ARG_NONE,
     NULL, OPT_CHECK_SPARSE,
     "check sparse matrix consistency", "check sparse"},
    {
     "outputfile", 'o', POPT_ARG_STRING,
     &outputfile_opt, OPT_OUTPUT_FILE,
     "output basename", "STRING"},
    {
     "importfile", 'i', POPT_ARG_STRING,
     &importfilename_opt, OPT_IMPORTFILE,
     "Import xml-data/(sparse,res,irm) files (comma separated list: -i f1,f2,...)",
     "FILE(S)"},
    {
     "logfile", 'l', POPT_ARG_STRING,
     &logfile_opt, OPT_LOG_FILE,
     "log filename", "FILE"},
    {
     "max-iter", '\0', POPT_ARG_INT,
     &max_iter_opt, OPT_MAX_ITER,
     "nb of maximum iteration ", "INT"},
    {
     "damping", 'd', POPT_ARG_FLOAT,
     &damping_opt, OPT_DAMPING,
     "damping factor", "FLOAT"},
    {
     "grad-damping", 'g', POPT_ARG_FLOAT,
     &grad_damping_opt, OPT_GRAD_DAMPING,
     "gradient damping factor", "FLOAT"},
    {
     "version", 'v', POPT_ARG_NONE,
     NULL, OPT_VERSION,
     "print version number", NULL},
    {NULL, '\0', 0, NULL, OPT_ZERO, NULL, NULL}
};

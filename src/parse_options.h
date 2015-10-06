#include <unistd.h>

#include <ray/ray3D.h>
#include <mesh/mesh.h>
#include <sparse/sparse.h>

#ifndef __LSQR_PARSE_H__
#define __LSQR_PARSE_H__

void parse_command_line(int argc, char **argv,
                        char **mesh_filename,
                        char **vmodel,
                        char **import_filename,
                        char **log_filename,
                        char **output_filename,
                        char **output_type,
                        int *max_iter,
                        double *damping,
                        double *grad_damping,
                        int *use_ach, int *check_sparse);

void lsqrsolve_info();

#endif

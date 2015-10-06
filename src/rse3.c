#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <signal.h>

#include <sparse/sparse.h>

#include "lsqr.h"
#include "lsqr_wrapper.h"
#include "catch_sig.h"
#include "extern.h"

/* rse3 : "resolution de systeme d'equations"
          this is just a test/example using sparse matrix (m x n) operator  + LSQR
	  sparse matrix file format is :
	
m n
i j val
...
*/

int please_stop_lsqr = 0;       /* stop lsqr and write the solution at the curent iteration */
int please_dump_lsqr = 0;       /* dump intermediate solution each  iterdump iteration */

/********/
/* MAIN */
/********/
int main(int argc, char *argv[])
{
    struct sparse_matrix_t *sparseA;
    struct vector_t *x, *b;

    lsqr_input *input;
    lsqr_output *output;
    lsqr_work *work;            /* zone temoraire de travail */
    lsqr_func *func;            /* func->mat_vec_prod -> APROD */

    /* cmd line arg */
    char *matrix_filename = NULL;
    char *vector_filename = NULL;
    char *sol_filename = NULL;
    int max_iter = -1;
    float damping = 0.;

    if (argc != 6 && argc != 7) {
        fprintf(stderr,
                "%s matrixfile vectorfile solutionfile nbitermax damping [iterdump]\n",
                argv[0]);
        exit(1);
    }
    matrix_filename = strdup(argv[1]);
    vector_filename = strdup(argv[2]);
    sol_filename = strdup(argv[3]);
    max_iter = (int) strtol(argv[4], (char **) NULL, 10);
    damping = strtod(argv[5], (char **) NULL);
    if (argc == 7) {
        please_dump_lsqr = (int) strtol(argv[6], (char **) NULL, 10);
    } else {
        please_dump_lsqr = 0;
    }

    /* read the sparse matrix */
    sparseA = read_ijk_sparse_matrix(matrix_filename, SPARSE_COL_LINK);
    fprintf(stderr, "read*matrix: ok (size=%ldx%ld, %ld elements)\n",
            sparseA->nb_line, sparseA->nb_col,
            sparseA->nb_line * sparseA->nb_col);
    show_sparse_stats(sparseA);

    /*write_sparse_matrix(sparseA, "test_sparse.txt"); */

    b = read_simple_vector(vector_filename);

    /*************************************************/
    /* check compatibility between matrix and vector */
    /*************************************************/
    if (sparseA->nb_line != b->length) {
        fprintf(stderr,
                "Error, check your matrix/vector sizes (%ld/%ld)\n",
                sparseA->nb_line, b->length);
        exit(1);
    }
    /* init vector solution to zero */
    x = new_vector(sparseA->nb_col);

    /* catch Ctrl-C signal */
    signal(SIGINT, emergency_halt);

    /*************************************************************/
    /* solve A.x = B                                             */
    /*************************************************************/

    /* LSQR alloc */
    alloc_lsqr_mem(&input, &output, &work, &func,
                   sparseA->nb_line, sparseA->nb_col);

    fprintf(stderr, "alloc_lsqr_mem : ok\n");

    /* defines the routine Mat.Vect to use */
    func->mat_vec_prod = sparseMATRIXxVECTOR;

    /* Set the input parameters for LSQR */
    input->num_rows = sparseA->nb_line;
    input->num_cols = sparseA->nb_col;
    input->rel_mat_err = .0;
    input->rel_rhs_err = .0;
    input->cond_lim = .0;
    input->lsqr_fp_out = stdout;
    input->rhs_vec = (dvec *) b;
    input->sol_vec = (dvec *) x;        /* initial guess */
    input->damp_val = damping;
    if (max_iter == -1) {
        input->max_iter = 4 * (sparseA->nb_col);
    } else {
        input->max_iter = max_iter;
    }

    /* resolution du systeme Ax=b */
    lsqr(input, output, work, func, sparseA);
    write_vector((struct vector_t *) output->sol_vec, sol_filename);

    /******************************************************/
    /* variance reduction, ie how the model fits the data */
    /* X = the final solution                             */
    /*                                                    */
    /*                 ||b-AX||²                          */
    /*         VR= 1 - --------                           */
    /*                  ||b||²                            */
    /*                                                    */
    /******************************************************/
    {
        double norm_b;
        double norm_b_AX;
        double VR;              /* variance reduction */

        struct vector_t *rhs;   /* right hand side */
        rhs = new_vector(sparseA->nb_line);

        /* use copy */
        dvec_copy((dvec *) b, (dvec *) rhs);

        norm_b = dvec_norm2((dvec *) rhs);

        /* does rhs = rhs + sparseA . output->sol_vec */
        /* here  rhs is overwritten */
        dvec_scale((-1.0), (dvec *) rhs);
        sparseMATRIXxVECTOR(0, output->sol_vec, (dvec *) rhs, sparseA);
        dvec_scale((-1.0), (dvec *) rhs);

        norm_b_AX = dvec_norm2((dvec *) rhs);

        //VR = 1 - (norm_b_AX*norm_b_AX)/(norm_b*norm_b);
        VR = 1 - norm_b_AX / norm_b;
        fprintf(stdout, "Variance reduction = %.2f%%\n", VR * 100);
        free_vector(rhs);
    }

    free_lsqr_mem(input, output, work, func);

    /* check A^t.A */
    /*
     * { struct sparse_matrix_t *AtA; AtA = AtransA (sparseA);
     * write_sparse_matrix(AtA, "AtA"); write_sparse_matrix(sparseA,
     * "A"); free_sparse_matrix (AtA);
     * 
     * } */

    free_sparse_matrix(sparseA);
    return (1);
}

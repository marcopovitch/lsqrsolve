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

/* rse : "resolution de systeme d'equations"
         this is just a test/example using sparse matrix operator (after sparsifying) + LSQR
*/

int please_stop_lsqr = 0;
int please_dump_lsqr = 0;

/********/

/* MAIN */

/********/
int main(int argc, char *argv[])
{
    struct matrix_t *matrixA;
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
    float damping = 0;

    if (argc != 4) {
        fprintf(stderr, "%s matrixfile vectorfile solutionfile\n",
                argv[0]);
        exit(1);
    }
    matrix_filename = strdup(argv[1]);
    vector_filename = strdup(argv[2]);
    sol_filename = strdup(argv[3]);

    /* read the  matrix */
    matrixA = read_matrix(matrix_filename);
    fprintf(stderr, "read*matrix: ok (size=%ldx%ld, %ld elements)\n",
            matrixA->nb_line, matrixA->nb_col,
            matrixA->nb_line * matrixA->nb_col);

    sparseA = sparsify(matrixA, SPARSE_COL_LINK);
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
    free_lsqr_mem(input, output, work, func);
    free_matrix(matrixA);

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

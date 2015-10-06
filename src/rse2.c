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

/* rse2 : "resolution de systeme d'equations"
         this is just a rse test/example using standard matrix operator + LSQR
*/

int please_stop_lsqr = 0;
int please_dump_lsqr = 0;

/********/

/* MAIN */

/********/
int main(int argc, char *argv[])
{
    struct matrix_t *matrixA;
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

    b = read_simple_vector(vector_filename);

        /*************************************************/
    /* check compatibility between matrix and vector */

        /*************************************************/
    if (matrixA->nb_line != b->length) {
        fprintf(stderr,
                "Error, check your matrix/vector sizes (%ld/%ld)\n",
                matrixA->nb_line, b->length);
        exit(1);
    }
    /* init vector solution to zero */
    x = new_vector(matrixA->nb_col);

    /* catch Ctrl-C signal */
    signal(SIGINT, emergency_halt);

        /*************************************************************/
    /* solve A.x = B                                             */

        /*************************************************************/

    /* LSQR alloc */
    alloc_lsqr_mem(&input, &output, &work, &func,
                   matrixA->nb_line, matrixA->nb_col);

    fprintf(stderr, "alloc_lsqr_mem : ok\n");

    /* defines the routine Mat.Vect to use */
    func->mat_vec_prod = MATRIXxVECTOR;

    /* Set the input parameters for LSQR */
    input->num_rows = matrixA->nb_line;
    input->num_cols = matrixA->nb_col;
    input->rel_mat_err = .0;
    input->rel_rhs_err = .0;
    input->cond_lim = .0;
    input->lsqr_fp_out = stdout;
    input->rhs_vec = (dvec *) b;
    input->sol_vec = (dvec *) x;        /* initial guess */
    input->damp_val = damping;
    if (max_iter == -1) {
        input->max_iter = 4 * (matrixA->nb_col);
    } else {
        input->max_iter = max_iter;
    }

    /* resolution du systeme Ax=b */
    lsqr(input, output, work, func, matrixA);
    write_vector((struct vector_t *) output->sol_vec, sol_filename);
    free_lsqr_mem(input, output, work, func);
    free_matrix(matrixA);
    return (1);
}

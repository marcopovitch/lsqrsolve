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

/* prod : compute A.X product   
          this is just a test/example using sparse matrix (m x n) operator
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
    struct vector_t *b;
    struct vector_t *rhs;       /* right hand side */

    /* cmd line arg */
    char *matrix_filename = NULL;
    char *vector_filename = NULL;
    char *sol_filename = NULL;

    if (argc != 4) {
        fprintf(stderr, "%s matrixfile vectorfile solutionfile\n",
                argv[0]);
        exit(1);
    }
    matrix_filename = strdup(argv[1]);
    vector_filename = strdup(argv[2]);
    sol_filename = strdup(argv[3]);

    /* read the sparse matrix */
    sparseA = read_ijk_sparse_matrix(matrix_filename, SPARSE_COL_LINK);
    fprintf(stderr, "read*matrix: ok (size=%ldx%ld, %ld elements)\n",
            sparseA->nb_line, sparseA->nb_col,
            sparseA->nb_line * sparseA->nb_col);
    show_sparse_stats(sparseA);

    b = read_simple_vector(vector_filename);

        /*************************************************/
    /* check compatibility between matrix and vector */
    /*  for product     Axb                          */

        /*************************************************/
    if (sparseA->nb_col != b->length) {
        fprintf(stderr,
                "Error, check your matrix/vector sizes (%ld/%ld)\n",
                sparseA->nb_col, b->length);
        exit(1);
    }

    fprintf(stderr, "Starting product ...\n");

    /* product */
    /*  re-use some code */
    /* compute if mode=0   rhs = rhs + A.b   */
    /* rhs is null */
    rhs = new_vector(sparseA->nb_line);
    sparseMATRIXxVECTOR(0, (dvec *) b, (dvec *) rhs, sparseA);

    write_vector((struct vector_t *) rhs, sol_filename);
    free_sparse_matrix(sparseA);
    /* free b, rhs */
    return (1);
}

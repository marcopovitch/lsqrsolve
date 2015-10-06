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

/* compute kernel density for Chris, using 
   sparse matrix file :
	
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
    struct vector_t *d;         /* density */
    int c;
    double accu;
    struct sparse_item_t *item;

    /* cmd line arg */
    char *matrix_filename = NULL;
    char *sol_filename = NULL;

    if (argc != 3) {
        fprintf(stderr, "%s matrixfile kernel_density_file\n", argv[0]);
        exit(1);
    }
    matrix_filename = strdup(argv[1]);
    sol_filename = strdup(argv[2]);

    /* read the sparse matrix */
    sparseA = read_ijk_sparse_matrix(matrix_filename, SPARSE_COL_LINK);
    fprintf(stderr, "read*matrix: ok (size=%ldx%ld, %ld elements)\n",
            sparseA->nb_line, sparseA->nb_col,
            sparseA->nb_line * sparseA->nb_col);
    show_sparse_stats(sparseA);

    /* allocation */
    d = new_vector(sparseA->nb_col);

    fprintf(stderr, "Starting kernel density computation ...\n");
    for (c = 0; c < sparseA->nb_col; c++) {
        item = sparseA->col[c];
        accu = 0.;
        while (item) {
            accu = accu + item->val;
            item = item->next_in_col;
        }
        d->mat[c] = accu;
    }
    write_vector((struct vector_t *) d, sol_filename);
    free_sparse_matrix(sparseA);
    /* free b, rhs */
    return (1);
}

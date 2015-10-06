#include "lsqr_wrapper.h"

const char *lsqr_msg[] = {
    "The exact solution is x = x0",
    "The residual Ax - b is small enough, given ATOL and BTOL",
    "The least squares error is small enough, given ATOL",
    "The estimated condition number has exceeded CONLIM",
    "The residual Ax - b is small enough, given machine precision",
    "The least squares error is small enough, given machine precision",
    "The estimated condition number has exceeded machine precision",
    "The iteration limit has been reached",
    "Interuption requested by the user (Ctrl-C)"
};

/*************************************/

/* computes :                        */

/* if mode = 0   ->  y = y + A.x   */

/* if mode = 1   ->  x = x + A^T.y */

/*************************************/
void MATRIXxVECTOR(long int mode, dvec * x, dvec * y, void *data)
{
    struct matrix_t *A;
    struct vector_t *tmp;
    long int i, j, k;

    A = (struct matrix_t *) data;

    if (mode == 0) {
        /* Compute  Y = Y + A*X */
        tmp = new_vector(A->nb_col);

        for (i = 0; i < A->nb_col; i++) {
            for (k = 0; k < A->nb_col; k++) {
                tmp->mat[i] += A->mat[i][k] * x->elements[k];
            }
        }

        for (j = 0; j < A->nb_col; j++) {
            y->elements[j] = y->elements[j] + tmp->mat[j];
        }
        return;
    }
    if (mode == 1) {
        /* Compute  X = X + A^T*Y */
        tmp = new_vector(A->nb_line);

        for (i = 0; i < A->nb_line; i++) {
            for (k = 0; k < A->nb_line; k++) {
                tmp->mat[i] += A->mat[k][i] * y->elements[k];
            }
        }

        for (j = 0; j < A->nb_line; j++) {
            x->elements[j] = x->elements[j] + tmp->mat[j];
        }

    }
}

/*************************************/

/* computes :                        */

/* if mode = 0   ->  y = y + A.x   */

/* if mode = 1   ->  x = x + A^T.y */

/*************************************/
void sparseMATRIXxVECTOR(long int mode, dvec * x, dvec * y, void *data)
{
    struct sparse_matrix_t *A;
    struct sparse_item_t *item;
    long int i, k;
    double accu;

    A = (struct sparse_matrix_t *) data;

    if (mode == 0) {
        /* Compute  Y = Y + A*X */
        for (i = 0; i < A->nb_line; i++) {
            item = A->line[i];
            accu = 0;
            while (item) {
                k = item->col_index;
                accu += item->val * x->elements[k];
                item = item->next_in_line;
            }
            y->elements[i] = y->elements[i] + accu;
        }

        return;
    }
    if (mode == 1) {
        /* Compute  X = X + A^T*Y */
        for (i = 0; i < A->nb_col; i++) {
            item = A->col[i];
            accu = 0.;
            while (item) {
                k = item->line_index;
                accu += item->val * y->elements[k];
                item = item->next_in_col;
            }
            x->elements[i] = x->elements[i] + accu;
        }
    }
}

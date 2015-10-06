#ifndef __REGUL_H__
#define __REGUL_H__

void regularization(struct mesh_t *mesh, char *sparsefile,
                    int nb_faces, double grad_damping,
                    long int start_line, long int *nb_lines);

void create_regul_DtD(struct sparse_matrix_t *A, long int *compress2fat,
                      struct mesh_t *mesh, char *sparsefile,
                      int nb_faces, double grad_damping,
                      long int *nb_lines);

void create_regul_DtD_irm(struct sparse_matrix_t *A,
                          long int *compress2fat, struct mesh_t *mesh,
                          char *sparsefile, int nb_faces,
                          double grad_damping, long int *nb_lines);

long int dicho_search_longint(long int *tab, long int ntab, long int key);

#endif

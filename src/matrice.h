#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifndef __MATRICE_H__
#define __MATRICE_H__

struct matrix_t {
    long int nb_line;
    long int nb_col;
    double **mat;
};

struct vector_t {
    long int length;
    double *mat;
};

struct matrix_t *new_matrix(int nb_line, int nb_col);
void free_matrix(struct matrix_t *c);

struct vector_t *new_vector(long int l);
void free_vector(struct vector_t *v);

void dump_matrix(char *txt, struct matrix_t *m);
void dump_vector(char *s, struct vector_t *v);

struct matrix_t *read_matrix(char *filename);
struct vector_t *read_vector(char *filename);
struct vector_t *read_simple_vector(char *filename);
struct vector_t *read_subvector(char *filename, long int *first,
                                long int *last);
struct vector_t *import_vector(struct vector_t *b, char *filename);
struct vector_t *vector_resize(struct vector_t *v, long int new_length);
void write_vector(struct vector_t *b, char *filename);

#endif

#include <stdio.h>
#include <stdlib.h>

#include <mesh/mesh.h>
#include <mesh/cell.h>

#include <sparse/sparse.h>

#ifndef  __IRM_H__
#define __IRM_H__

struct irm_t {
    int x;
    int y;
    int z;
    int m;
    int n;
    float s;
};

struct irm_t *read_irm(char *filename, int *nb_metacell);
void irm_update(struct vector_t *solution, struct irm_t **irm,
                int *nb_metacell, int nb_irm, struct mesh_t *mesh);
void free_irm(struct irm_t **irm, int nb_irm);

#endif

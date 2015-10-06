#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <locale.h>

#include <mesh/mesh.h>
#include <mesh/layer.h>
#include <mesh/cell.h>
#include <mesh/cellinfo.h>

#include <ray/vitesse.h>

#include <sparse/sparse.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef __EXPORT_H__
#define __EXPORT_H__

void export2matlab(struct vector_t *v, char *filename,
                   struct mesh_t *mesh, struct velocity_model_t *vm,
                   long int nb_iter,
                   double damping, double grad_damping, int use_ach);

void export2gmt(struct vector_t *v, char *filename,
                struct mesh_t *mesh, struct velocity_model_t *vm,
                long int nb_iter,
                double damping, double grad_damping, int use_ach);

void export2sco(struct vector_t *v, char *filename, struct mesh_t *mesh,
                struct velocity_model_t *vm, long int nb_iter,
                double damping, double grad_damping, int use_ach);

#endif

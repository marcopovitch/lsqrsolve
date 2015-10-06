#include "export.h"

/* the %velocity is computed in respect to 
 * the center of the layer 
 */

#define VREF "mean"

void export2matlab(struct vector_t *v,
                   char *filename,
                   struct mesh_t *mesh,
                   struct velocity_model_t *vm,
                   long int nb_iter,
                   double damping, double grad_damping, int use_ach)
{
    long int i;
    int layer;
    FILE *fd;
    int nbx, nby;
    char *matlab_file;
    int x, y, z;
    double vel, delta_v, delta_s, delta_v_percent;
    char *field[3] = { "slowness", "%velocity" };
    char *pfield;

    if (use_ach) {
        pfield = field[1];
    } else if (!vm) {
        pfield = field[0];
    } else {
        pfield = field[1];
    }

    matlab_file = (char *)
        malloc((strlen(filename) + strlen(".m") + 1) * sizeof(char));
    assert(matlab_file);

    sprintf(matlab_file, "%s.m", filename);
    fprintf(stdout, "Writing layer ouput to matlab file '%s'\n",
            matlab_file);
    fprintf(stdout, "data type is %s, %ld items\n", pfield, v->length);

    if (!(fd = fopen(matlab_file, "w"))) {
        perror(matlab_file);
        exit(1);
    }

    layer = 0;
    for (i = 0; i < v->length; i++) {

        unlinearize_cell_id(i, &x, &y, &z, mesh);
        if (i == 0 || layer != z) {
            /* new layer : write new header */
            layer = z;
            if (layer != 0) {
                fprintf(fd, "\n");
            }

            nbx = mesh->layer[layer]->nlat;
            nby = mesh->layer[layer]->nlon;

            fprintf(fd,
                    "# %s v%s, mesh=%s, nb_iter=%ld, damping(d/g)=%f/%f, data=%s, vref=%s\n",
                    PACKAGE, VERSION, mesh->xml_filename, nb_iter, damping,
                    grad_damping, pfield, VREF);
            fprintf(fd, "# name: layer%d\n", layer);
            fprintf(fd, "# type: matrix\n");
            fprintf(fd, "# rows: %d\n", nbx);
            fprintf(fd, "# columns: %d\n", nby);
        }

        if (use_ach) {
            fprintf(fd, "%g ", v->mat[i] * -100.0);
        } else if (!vm) {
            /* slowness */
            fprintf(fd, "%g ", v->mat[i]);
        } else {
            vel = compute_vmean(vm, mesh->layer[z]->zstart,
                                mesh->layer[z]->zend);
            delta_s = v->mat[i];
            delta_v = 1 / (1 / vel + delta_s) - vel;
            delta_v_percent = delta_v / vel * 100;
            fprintf(fd, "%g ", delta_v_percent);
        }
    }
    fclose(fd);
    free(matlab_file);
}

void export2gmt(struct vector_t *v,
                char *filename,
                struct mesh_t *mesh,
                struct velocity_model_t *vm,
                long int nb_iter,
                double damping, double grad_damping, int use_ach)
{
    long int i;
    int layer;
    FILE *fd = NULL;
    char *layerfile = NULL;
    int nbx, nby;
    int x, y, z;
    double vel, delta_v, delta_s, delta_v_percent;
    char *field[3] =
        { "slowness", "%velocity",
"slowness & %velocity & mean_velocity" };
    char *pfield;

    if (use_ach) {
        pfield = field[1];
    } else if (!vm) {
        pfield = field[0];
    } else {
        pfield = field[2];
    }

    fprintf(stdout, "Writing layer ouput to GMT file '%s-*' ", filename);
    fprintf(stdout, "data type is %s, %ld items\n", pfield, v->length);

    layer = 0;
    nbx = 0;
    nby = 0;

    for (i = 0; i < v->length; i++) {

        unlinearize_cell_id(i, &x, &y, &z, mesh);

        if (i == 0 || layer != z) {
            /* start a new layer */
            layer = z;
            nbx = mesh->layer[layer]->nlat;
            nby = mesh->layer[layer]->nlon;

            if (fd)
                fclose(fd);
            layerfile = (char *)
                realloc(layerfile,
                        (strlen(filename) + strlen(".xy") +
                         4) * sizeof(char));
            assert(layerfile);
            sprintf(layerfile, "%s-%.2d.xy", filename, layer);
            fprintf(stdout, "Writing to file %s\n", layerfile);

            if (!(fd = fopen(layerfile, "w"))) {
                perror(layerfile);
                exit(1);
            }
            fprintf(fd,
                    "# %s v%s, mesh=%s, nb_iter=%ld, damping(d/g)=%f/%f\n",
                    PACKAGE, VERSION, mesh->xml_filename, nb_iter, damping,
                    grad_damping);
            fprintf(fd, "# %s : data is %s, layer %d, nlat=%d, nlon=%d\n",
                    layerfile, pfield, layer, nbx, nby);
        }

        if (use_ach) {
            fprintf(fd, "%g\n", v->mat[i] * -100.0);
        } else if (!vm) {
            fprintf(fd, "%g\n", v->mat[i]);
        } else {
            vel = compute_vmean(vm, mesh->layer[z]->zstart,
                                mesh->layer[z]->zend);
            delta_s = v->mat[i];
            delta_v = 1 / (1 / vel + delta_s) - vel;
            delta_v_percent = delta_v / vel * 100;
            fprintf(fd, "%10e %10e %f\n", delta_s, delta_v_percent, vel);
        }
    }
    free(layerfile);
    fclose(fd);
}

void export2sco(struct vector_t *v,
                char *filename,
                struct mesh_t *mesh,
                struct velocity_model_t *vm,
                long int nb_iter,
                double damping, double grad_damping, int use_ach)
{
    long int i;
    FILE *fd;
    int x, y, z = -1;
    int nb_score;
    char *sco_file;

    /* vel */
    double vel, delta_v, delta_s, delta_v_percent;

    sco_file = (char *)
        malloc((strlen(filename) + strlen(".sco") + 1) * sizeof(char));
    assert(sco_file);

    sprintf(sco_file, "%s.sco", filename);
    fprintf(stdout, "Writing layer ouput to sco file '%s'\n", sco_file);

    if (!(fd = fopen(sco_file, "w"))) {
        perror(sco_file);
        exit(1);
    }

    fprintf(fd,
            "# format=sco, generated by %s v%s, mesh=%s, nb_iter=%ld, damping(d/g)=%f/%f\n",
            PACKAGE, VERSION, mesh->xml_filename, nb_iter, damping,
            grad_damping);
    if (use_ach) {
        nb_score = 1;
        fprintf(fd, "%d delta_velocity\n", nb_score);
    } else if (!vm) {
        nb_score = 1;
        fprintf(fd, "%d delta_slowness \n", nb_score);
    } else {
        nb_score = 3;
        fprintf(fd, "%d delta_slowness delta_velocity mean_velocity\n",
                nb_score);
    }

    /* export data */
    for (i = 0; i < v->length; i++) {
        unlinearize_cell_id(i, &x, &y, &z, mesh);
        if (fabs(v->mat[i]) > EPS_SPARSE) {
            if (use_ach) {
                fprintf(fd, "[%d,%d,%d] %g\n", x, y, z,
                        v->mat[i] * -100.0);
            } else if (!vm) {
                fprintf(fd, "[%d,%d,%d] %g\n", x, y, z, v->mat[i]);
            } else {
                vel = compute_vmean(vm, mesh->layer[z]->zstart,
                                    mesh->layer[z]->zend);
                delta_s = v->mat[i];
                delta_v = 1 / (1 / vel + delta_s) - vel;
                delta_v_percent = delta_v / vel * 100;
                fprintf(fd, "[%d,%d,%d] %10e %10e %f\n",
                        x, y, z, delta_s, delta_v_percent, vel);
            }
        }
    }
    fclose(fd);

    /* update the xml data structure */
    mesh_remove_data_entry(mesh, SCO);
    mesh_add_data_filename(mesh, SCO, sco_file);

    free(sco_file);
}

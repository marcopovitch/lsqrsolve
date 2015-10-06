#include "irm.h"

/** \brief read a irm file and return an irm_t array null terminated */
struct irm_t *read_irm(char *filename, int *nb_metacell)
{
    FILE *fd;
    struct irm_t *irm;
    int nb_read, i;
    int cpt = 0;
    char line[256];

    fprintf(stdout, "reading irregular mesh  from '%s' ... ", filename);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }

    /* skip 3 lines of header */
    fgets(line, 256, fd);
    fgets(line, 256, fd);
    fgets(line, 256, fd);

    nb_read = fscanf(fd, "%d\n", nb_metacell);
    if (nb_read != 1) {
        fprintf(stdout, "\nError reading nb_metacell in '%s'\n", filename);
        exit(1);
    }

    irm = (struct irm_t *) malloc(sizeof(struct irm_t) * (*nb_metacell));
    assert(irm);

    for (i = 0; i < *nb_metacell; i++) {
        if (feof(fd)) {
            break;
        }

        nb_read = fscanf(fd, "%d %d %d %d %d %f\n",
                         &(irm[i].x), &(irm[i].y), &(irm[i].z),
                         &(irm[i].m), &(irm[i].n), &(irm[i].s));

        if (nb_read != 6) {
            fprintf(stdout, "\n");
            fprintf(stderr, "read_irm: file '%s' corrupted\n", filename);
            exit(1);
        }

        cpt++;
    }

    if (cpt != *nb_metacell) {
        fprintf(stdout, "\n");
        fprintf(stderr, "read_irm: file '%s' corrupted read %d/%d\n",
                filename, cpt, *nb_metacell);
        exit(1);
    }

    fprintf(stdout, "loaded %d metacells\n", *nb_metacell);
    fflush(stdout);

    fclose(fd);
    return (irm);

}

void irm_update(struct vector_t *solution, struct irm_t **irm,
                int *nb_metacell, int nb_irm, struct mesh_t *mesh)
{
    struct coord_z3_t cell, orig_cell;
    int i, mc;
    long int offset;
    int mm, nn;
    int m, n;
    double val;
    struct irm_t *cur_irm;

    /* x=lat y=lon z=layer */
    /* m follows y */
    /* n follows x */

    fprintf(stdout, "Update solution using irm  ... ");

    for (i = 0; i < nb_irm; i++) {
        cur_irm = irm[i];
        for (mc = 0; mc < nb_metacell[i]; mc++) {
            orig_cell.x = cell.x = cur_irm[mc].x;
            orig_cell.y = cell.y = cur_irm[mc].y;
            orig_cell.z = cell.z = cur_irm[mc].z;
            m = irm[i][mc].m;
            n = irm[i][mc].n;
            offset = linearize_cell_id(&cell, mesh);
            val = solution->mat[offset];

            /* update the solution */
            for (mm = 0; mm < m; mm++) {
                cell.x = orig_cell.x;
                cell.y = orig_cell.y + mm;
                cell.z = orig_cell.z;

                for (nn = 0; nn < n; nn++) {
                    offset = linearize_cell_id(&cell, mesh);
                    solution->mat[offset] = val;
                    cell.x--;
                }
            }
        }
    }

    fprintf(stdout, "ok\n");
}

void free_irm(struct irm_t **irm, int nb_irm)
{
    int i;

    for (i = 0; i < nb_irm; i++) {
        free(irm[i]);
    }
    free(irm);

}

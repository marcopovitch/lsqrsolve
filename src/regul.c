#include <mesh/mesh.h>
#include <mesh/layer.h>
#include <mesh/cell.h>
#include <mesh/metacell.h>

#include <sparse/sparse.h>

#include "regul.h"

/** \brief create a sparse matrix to regularize a linear system using the
 *   gradient damping.
 *  
 **/
void regularization(struct mesh_t *mesh, char *sparsefile,
                    int nb_faces, double grad_damping,
                    long int start_line, long int *nb_lines)
{
    int nlayer;
    struct sparse_matrix_t *regul, *AtA;
    struct cell_t *current_cell, *last_cell = NULL;
    struct cell_t *first_cell_on_parallel, *neighbour;
    int ncell = 0, nparallel = 0;
    int north;
    long int lcid, neighbour_lcid, current_line;
    int face;
    int face_list[2] = { 0 /*NSD_1 */ , 2 /*EAST*/ };

    if (!mesh) {
        fprintf(stderr, "regularization: can't initialize mesh\n");
        exit(1);
    }

    make_mesh(mesh);

    regul =
        new_sparse_matrix(2 * mesh->ncells, mesh->ncells, SPARSE_COL_LINK);
    fprintf(stderr, "regularization new_sparse_matrix (%p) (%ldx%ld)\n",
            regul, regul->nb_line, regul->nb_col);

    current_line = 0;
    for (nlayer = 0; nlayer < mesh->nlayers; nlayer++) {

        north = get_north_direction(NULL, mesh->layer[nlayer], mesh);
        current_cell = mesh->layer[nlayer]->cell;

        do {
            nparallel++;
            first_cell_on_parallel = current_cell;

            /* go to the EAST as long as we can, but don't loop forever */
            do {

                lcid = linearize_cell_id(&(current_cell->id), mesh);

                /* get neighbour lcid */
                for (face = 0; face < 2; face++) {
                    neighbour =
                        get_cell_from_list(current_cell->
                                           neighbour_list[face_list[face]],
                                           0);
                    if (!neighbour) {
                        continue;
                    }

                    neighbour_lcid =
                        linearize_cell_id(&(neighbour->id), mesh);
                    sparse_set_value(regul, current_line, neighbour_lcid,
                                     -1. * grad_damping, NULL);
                    sparse_set_value(regul, current_line, lcid,
                                     grad_damping, NULL);
                    current_line++;
                }

                /* next cell */
                last_cell = current_cell;
                current_cell =
                    get_cell_from_list(current_cell->
                                       neighbour_list[EAST_D], 0);
                ncell++;

            } while (current_cell
                     && current_cell != first_cell_on_parallel);

            /* go to the north, as long as we can, ie. next parallel */
            current_cell =
                get_cell_from_list(first_cell_on_parallel->
                                   neighbour_list[north], 0);

        } while (current_cell
                 && !cells_are_in_the_same_crown(current_cell, last_cell));
    }

    *nb_lines = current_line;

    /*fprintf(stderr, "regularization: %d cells on %d layers\n", ncell, nlayer); */

    AtA = AtransA(regul);

    /* show output */
    /*write_sparse_matrix(AtA, "DtD.sparse");
       write_sparse_matrix(regul, "D.sparse"); */
    free_sparse_matrix(regul);

    write_sparse_matrix_with_line_offset(AtA, start_line, sparsefile);
    free_sparse_matrix(AtA);
}

/** \brief constructs the regularization matrix DtD based on smooth
 *  slowness gradient between neighbours cells
 *
 * When USE_COMPRESSION is activated, 
 * all empties cells are virtually removed
 * from the inversion. Those cells will not
 * have any associated velocity in the output files.
 * This will also modified the regularization
 * behaviour, such as gradient smoothing is only done
 * with non empty cells.
 *
 *  DtD.s = 0
 */
void create_regul_DtD(struct sparse_matrix_t *A, long int *compress2fat,
                      struct mesh_t *mesh, char *sparsefile,
                      int nb_faces, double grad_damping,
                      long int *nb_lines)
{
    int nlayer;
    struct sparse_matrix_t *regul;
    struct cell_t *current_cell, *last_cell = NULL;

    struct cell_t *first_cell_on_parallel, *neighbour;
    int ncell = 0, nparallel = 0, nb_neighbours = 0;
    int north;
    long int lcid, neighbour_lcid, current_line;
    int face;

    if (!mesh) {
        fprintf(stderr, "%s: can't initialize mesh\n", __FUNCTION__);
        exit(1);
    }

    make_mesh(mesh);

    regul = new_sparse_matrix(A->nb_col, A->nb_col, SPARSE_COL_LINK);
    fprintf(stderr, "%s: new_sparse_matrix (%p) (%ldx%ld)\n",
            __FUNCTION__, regul, regul->nb_line, regul->nb_col);
    fprintf(stderr, "%s: using %d faces\n", __FUNCTION__, nb_faces);

    current_line = 0;
    for (nlayer = 0; nlayer < mesh->nlayers; nlayer++) {

        north = get_north_direction(NULL, mesh->layer[nlayer], mesh);
        current_cell = mesh->layer[nlayer]->cell;

        do {
            nparallel++;
            first_cell_on_parallel = current_cell;

            /* go to the EAST as long as we can, but don't loop forever */
            do {

                lcid = linearize_cell_id(&(current_cell->id), mesh);

                /* check if current_cell carry ray information 
                 * if not, no regularization done.
                 */
                if (compress2fat) {
                    lcid =
                        dicho_search_longint(compress2fat, A->nb_col,
                                             lcid);
                    if (lcid == -1 || !A->col[lcid]) {
                        goto NEXT_CELL;
                    }
                }

                /* get neighbour lcid */
                nb_neighbours = 0;
                for (face = 0; face < nb_faces; face++) {
                    neighbour =
                        get_cell_from_list(current_cell->
                                           neighbour_list[face], 0);

                    if (!neighbour) {
                        continue;
                    }

                    neighbour_lcid =
                        linearize_cell_id(&(neighbour->id), mesh);

                    /* same as before: check if current_cell carry ray information 
                     * if not, no regularization done.
                     */
                    if (compress2fat) {
                        neighbour_lcid =
                            dicho_search_longint(compress2fat, A->nb_col,
                                                 neighbour_lcid);
                        if (neighbour_lcid == -1
                            || !A->col[neighbour_lcid]) {
                            continue;
                        }
                    }

                    /* do the regularization */
                    sparse_set_value(regul, current_line, neighbour_lcid,
                                     -1. * grad_damping, NULL);
                    nb_neighbours++;
                }

                sparse_set_value(regul, current_line, lcid,
                                 nb_neighbours * grad_damping, NULL);
                current_line++;

                /* next cell */
              NEXT_CELL:last_cell = current_cell;
                current_cell =
                    get_cell_from_list(current_cell->
                                       neighbour_list[EAST_D], 0);
                ncell++;

            } while (current_cell
                     && current_cell != first_cell_on_parallel);

            /* go to the north, as long as we can, ie. next parallel */
            current_cell =
                get_cell_from_list(first_cell_on_parallel->
                                   neighbour_list[north], 0);

        } while (current_cell
                 && !cells_are_in_the_same_crown(current_cell, last_cell));
    }

    *nb_lines = current_line;

    /*write_sparse_matrix(regul, "DtD2.sparse"); */
    write_sparse_matrix_with_line_offset(regul, A->nb_line, sparsefile);

    free_sparse_matrix(regul);
}

void create_regul_DtD_irm(struct sparse_matrix_t *A,
                          long int *compress2fat, struct mesh_t *mesh,
                          char *sparsefile, int nb_faces,
                          double grad_damping, long int *nb_lines)
{
    int nlayer, i;
    int nb_neighbours, nb_total_neighbours;
    int face, neighb;
    long int old_lcid, old_ngh_lcid, lcid, neighbour_lcid, current_line;
    struct sparse_matrix_t *regul;
    struct cell_t *head_cell, *neighbour_cell;

    if (!mesh) {
        fprintf(stderr, "%s: can't initialize mesh\n", __FUNCTION__);
        exit(1);
    }

    if (!mesh->nb_total_metacell) {
        fprintf(stderr, "%s: irm files are not yet loaded !\n",
                __FUNCTION__);
        exit(1);
    }

    regul = new_sparse_matrix(A->nb_col, A->nb_col, SPARSE_COL_LINK);
    fprintf(stderr, "%s: new_sparse_matrix (%p) (%ldx%ld)\n",
            __FUNCTION__, regul, regul->nb_line, regul->nb_col);
    /*fprintf(stderr, "%s: using %d faces\n", __FUNCTION__, nb_faces); */

    current_line = 0;
    for (nlayer = 0; nlayer < mesh->nlayers; nlayer++) {
        fprintf(stderr, "%s: working on layer %d\n", __FUNCTION__, nlayer);
        for (i = 0; i < mesh->nb_metacell[nlayer]; i++) {

            /* the heading cell */
            head_cell = mesh->metacell[nlayer][i];
            old_lcid = lcid = linearize_cell_id(&(head_cell->id), mesh);
            lcid = dicho_search_longint(compress2fat, A->nb_col, lcid);

            if (lcid == -1) {
                fprintf(stderr,
                        "%s: can't find Hcellid=%ld in compressed matrix (nbcol=%ld)\n",
                        __FUNCTION__, old_lcid, A->nb_col);
                fprintf(stderr,
                        "%s: FIXME, only compress column if:\n\t1. no ray (done), 2. no head cell (to be done).\n",
                        __FUNCTION__);
                exit(1);
            }

            /* neighbours on the 4 faces N/S/E/W */
            nb_total_neighbours = 0;
            for (face = 0; face < 4; face++) {
                nb_neighbours =
                    cell_list_get_nb_item(head_cell->
                                          meta_neighbour_list[face]);

                for (neighb = 0; neighb < nb_neighbours; neighb++) {
                    neighbour_cell =
                        get_cell_from_list(head_cell->
                                           meta_neighbour_list[face],
                                           neighb);
                    old_ngh_lcid = neighbour_lcid =
                        linearize_cell_id(&(neighbour_cell->id), mesh);
                    neighbour_lcid =
                        dicho_search_longint(compress2fat, A->nb_col,
                                             neighbour_lcid);

                    if (neighbour_lcid == -1) {
                        fprintf(stderr,
                                "%s: can't find neighbour_cellid=%ld in compressed matrix (nbcol=%ld)\n",
                                __FUNCTION__, old_ngh_lcid, A->nb_col);
                        fprintf(stderr,
                                "%s: neighbour_cell=%p score=%.2f\n",
                                __FUNCTION__, neighbour_cell,
                                neighbour_cell->meta_scores[0].score);
                        /* the neighbour cell is not crossed by any ray: skip it */
                    } else {

                        /* do the regularization */
                        nb_total_neighbours++;
                        sparse_set_value(regul, current_line,
                                         neighbour_lcid,
                                         -1. * grad_damping, NULL);

                        if (old_lcid == old_ngh_lcid) {
                            fprintf(stderr,
                                    "old_lcid=%ld == old_ngh_lcid=%ld\n",
                                    old_lcid, old_ngh_lcid);
                        }
                    }

                }
            }

            sparse_set_value(regul, current_line, lcid,
                             nb_total_neighbours * grad_damping, NULL);
            current_line++;
        }
    }

    *nb_lines = current_line;

    write_sparse_matrix(regul, "DtD2-irm.sparse");
    write_sparse_matrix_with_line_offset(regul, A->nb_line, sparsefile);

    free_sparse_matrix(regul);
}

long int dicho_search_longint(long int *tab, long int ntab, long int key)
{
    long int low, high;
    long int i;

    assert(tab);

    for (low = -1, high = ntab; high - low > 1;) {
        i = (high + low) / 2;

        if (key <= tab[i]) {
            high = i;
        } else {
            low = i;
        }
    }
    if (key == tab[high]) {
        return (high);
    } else {
        return (-1);
    }
}

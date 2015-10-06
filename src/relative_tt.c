#include "relative_tt.h"

/* \brief import seismic event file */
struct event_list_t *import_event_file(char *filename)
{
    FILE *fd;
    int nb_read;
    long int event_id, ray_id, last_event_id, nb_rays, i;
    double time;
    struct event_list_t *evl;
    struct event_info_t *event_info = NULL;

    fprintf(stdout, "reading event file from '%s' ... ", filename);
    fflush(stdout);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }

    /* event list alloc */
    evl = (struct event_list_t *) calloc(1, sizeof(struct event_list_t));
    assert(evl);

    last_event_id = -1;
    while (!feof(fd)) {

        nb_read = fscanf(fd, "%ld %ld %lf\n", &event_id, &ray_id, &time);
        if (nb_read != 3) {
            fprintf(stdout, "\n");
            fprintf(stderr, "Error reading in '%s'\n", filename);
            exit(1);
        }

        /* store */
        if (event_id != last_event_id) {
            /* new event */
            event_info = (struct event_info_t *)
                calloc(1, sizeof(struct event_info_t));
            assert(event_info);
            event_info->event_id = event_id;

            /* update the list */
            evl->nb_event++;
            evl->event_array = (struct event_info_t **)
                realloc(evl->event_array,
                        sizeof(struct event_info_t *) * evl->nb_event);
            assert(evl->event_array);
            evl->event_array[evl->nb_event - 1] = event_info;
        }

        /* add item */
        event_info->nb_ray++;
        event_info->raytime = (struct rayid_time_t *)
            realloc(event_info->raytime,
                    event_info->nb_ray * sizeof(struct rayid_time_t));
        assert(event_info->raytime);
        event_info->raytime[event_info->nb_ray - 1].rayid = ray_id;
        event_info->raytime[event_info->nb_ray - 1].time = time;

        last_event_id = event_id;
    }

    nb_rays = 0;
    for (i = 0; i < evl->nb_event; i++) {
        nb_rays += evl->event_array[i]->nb_ray;
    }

    fprintf(stdout, "loaded %ld events, %ld rays\n", evl->nb_event,
            nb_rays);

    fflush(stdout);

    return (evl);
}

/* \brief ACH teleseismic tomographic inversion */
void relative_tt(struct sparse_matrix_t *M, struct vector_t *b,
                 char *filename)
{
    long int i, r, n, rayid;
    struct sparse_item_t *cur_item;
    struct event_list_t *evl;
    struct event_info_t *event_info;
    double *accu;

    evl = import_event_file(filename);
    fprintf(stderr, "using relative traveltime ... working on event ");

    /* alloc */
    accu = (double *) calloc(M->nb_col, sizeof(double));
    assert(accu);

    /* for a given event 'i' compute the average residual 
       over all station that have recorded this event, 
       then update b vector accordingly               */
    for (i = 0; i < evl->nb_event; i++) {
        event_info = evl->event_array[i];
        fprintf(stderr, "[%-8ld]\b\b\b\b\b\b\b\b\b\b",
                event_info->event_id);
        /*fprintf(stderr, "event=%ld nbray=%ld\n",  
           event_info->event_id, event_info->nb_ray); */

        /* compute residual mean */
        event_info->mean = 0.;
        for (r = 0; r < event_info->nb_ray; r++) {
            rayid = event_info->raytime[r].rayid;
            event_info->mean += b->mat[rayid];
        }
        event_info->mean = event_info->mean / event_info->nb_ray;
        /*fprintf(stderr, "mean[%ld]=%f\n", event_info->event_id, event_info->mean); */

        /* update b */
        for (r = 0; r < event_info->nb_ray; r++) {
            rayid = event_info->raytime[r].rayid;
            b->mat[rayid] -= event_info->mean;
            /*fprintf(stderr, "b[%ld]=%f\n", rayid, b->mat[rayid]); */
        }

        /* for a given event 'i', 
         * remove for each block and  each ray-paths (which linked all the stations to event 'i')
         * the mean of those ray-path concerning event 'i' in the same block */

        /* reset accu */
        for (n = 0; n < M->nb_col; n++) {
            accu[n] = 0.;
        }

        /* update M */
        for (r = 0; r < event_info->nb_ray; r++) {
            rayid = event_info->raytime[r].rayid;

            cur_item = M->line[rayid];
            while (cur_item) {
                n = cur_item->col_index;
                accu[n] += cur_item->val;
                cur_item = cur_item->next_in_line;
            }
        }

        for (r = 0; r < event_info->nb_ray; r++) {
            rayid = event_info->raytime[r].rayid;
            /*fprintf(stderr, "rayid=%ld\n", rayid); */

            cur_item = M->line[rayid];
            while (cur_item) {
                n = cur_item->col_index;
                /*fprintf(stderr, "\tM[%ld,%ld]=%f  ", 
                   rayid, n, cur_item->val); */
                cur_item->val -= accu[n] / event_info->nb_ray;

                /*fprintf(stderr, "M[%ld,%ld]=%f accu[%ld]=%f\n", 
                   rayid, n, cur_item->val, 
                   n, accu[n]); */

                cur_item = cur_item->next_in_line;
            }
        }

    }

    fprintf(stderr, "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
    fprintf(stderr, "%ld events processed        \n", evl->nb_event);

    /* free */
    free(accu);
    for (i = 0; i < evl->nb_event; i++) {
        event_info = evl->event_array[i];
        free(event_info->raytime);
        free(event_info);
    }
    free(evl);
}

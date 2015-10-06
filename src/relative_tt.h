#include <sparse/sparse.h>

#ifndef __REALATIVE_TT_H__
#define __REALATIVE_TT_H__
struct rayid_time_t {
    long int rayid;
    double time;
};

struct event_info_t {
    long int event_id;
    long int nb_ray;
    double mean;
    struct rayid_time_t *raytime;
};

struct event_list_t {
    long int nb_event;
    struct event_info_t **event_array;
};

struct event_list_t *import_event_file(char *filename);
void relative_tt(struct sparse_matrix_t *M, struct vector_t *b,
                 char *filename);
#endif

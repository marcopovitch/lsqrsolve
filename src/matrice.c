#include "matrice.h"

/*********************************/

/* create a new matrix with zero */

/*********************************/
struct matrix_t *new_matrix(int nb_line, int nb_col)
{
    struct matrix_t *c;
    long int i;

    c = (struct matrix_t *) malloc(sizeof(struct matrix_t));
    assert(c);

    c->nb_line = nb_line;
    c->nb_col = nb_col;
    c->mat = (double **) malloc(c->nb_line * sizeof(double *));

    for (i = 0; i < c->nb_line; i++) {
        (c->mat)[i] = (double *) calloc(c->nb_col, sizeof(double));
        assert((c->mat)[i]);
    }

    return (c);
}

/***************/

/* free matrix */

/***************/
void free_matrix(struct matrix_t *c)
{
    long int i;
    for (i = 0; i < c->nb_line; i++) {
        free((c->mat)[i]);
    }
    free(c->mat);
    free(c);
    c = NULL;
}

void dump_matrix(char *txt, struct matrix_t *m)
{
    long int i, j;

    if (m == NULL) {
        fprintf(stderr, "%s : matrix is not defined\n", txt);
        return;
    }
    fprintf(stderr, "matrix %s : addr=%p  (%ld x %ld)\n", txt, m,
            m->nb_line, m->nb_col);
    for (i = 0; i < m->nb_line; i++) {
        fprintf(stderr, "\t");
        for (j = 0; j < m->nb_col; j++) {
            fprintf(stderr, "%.2f ", (m->mat)[i][j]);
        }
        fprintf(stderr, "\n");
    }
}

/*********************************/

/* create a new vector with zero */

/*********************************/
struct vector_t *new_vector(long int l)
{
    struct vector_t *v;

    v = (struct vector_t *) malloc(sizeof(struct vector_t));
    assert(v);
    v->mat = (double *) calloc(l, sizeof(double));
    assert(v->mat);
    v->length = l;
    return (v);
}

void free_vector(struct vector_t *v)
{
    free(v->mat);
    free(v);
    v = NULL;
}

void dump_vector(char *s, struct vector_t *v)
{
    long int i;

    fprintf(stderr, "vector %s : ", s);
    for (i = 0; i < v->length; i++)
        fprintf(stderr, "%f ", v->mat[i]);
    fprintf(stderr, "\n");
}

/********************************************************/

/* Chargement de matrice et vecteur a partir de fichier */

/********************************************************/
struct matrix_t *read_matrix(char *filename)
{
    struct matrix_t *a;
    int m, n, i, j;
    FILE *fd;
    int nb_read;

    fprintf(stdout, "reading matrix A from file '%s'\n", filename);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }

    nb_read = fscanf(fd, "%d", &m);
    if (nb_read != 1) {
        fprintf(stderr, "Error reading 'm' in '%s'\n", filename);
        exit(1);
    }

    nb_read = fscanf(fd, "%d", &n);
    if (nb_read != 1) {
        fprintf(stderr, "Error reading 'n' in '%s'\n", filename);
        exit(1);
    }

    a = new_matrix(m, n);

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            nb_read = fscanf(fd, "%lf", &(a->mat[i][j]));
            if (nb_read != 1) {
                fprintf(stderr, "Error reading mat[%d][%d] in '%s'\n",
                        i, j, filename);
                exit(1);
            }
        }
    }

    fclose(fd);

    return (a);
}

struct vector_t *import_vector(struct vector_t *b, char *filename)
{
    long int n, j;
    FILE *fd;
    int nb_read;
    long int rayid;
    double val;

    if (!b) {
        b = read_vector(filename);
        return (b);
    }

    fprintf(stdout, "importing vector b from '%s' into (%p) ... ",
            filename, b);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }

    nb_read = fscanf(fd, "%ld", &n);
    if (nb_read != 1) {
        fprintf(stderr, "Error reading 'n' in '%s'\n", filename);
        exit(1);
    }

    if (n != b->length) {
        fprintf(stderr,
                "import_vector failed: size required %ld, read %ld in %s\n",
                b->length, n, filename);
        exit(1);
    }

    /* load vector data from file */
    j = 0;
    while (!feof(fd)) {
        nb_read = fscanf(fd, "%ld %lf\n", &rayid, &val);

        if (nb_read != 2) {
            fprintf(stderr, "Error reading mat[%ld] in '%s'\n",
                    j, filename);
            exit(1);
        }

        if (fabs(b->mat[rayid]) > 1.0e-6) {
            fprintf(stdout,
                    "import_vector: duplicate value (%ld) old=%f/new=%f\n",
                    rayid, b->mat[rayid], val);
        }
        b->mat[rayid] = val;
        j++;
    }
    fclose(fd);
    fprintf(stdout, "%ld lines\n", j);
    fflush(stdout);
    return (b);

}

/** \brief read a file as vector
 * 
 * the file is formated as follow :
 *
 * nbitem
 * value1
 * value2
 * ...
 */
struct vector_t *read_simple_vector(char *filename)
{
    struct vector_t *b;
    long int n, j;
    FILE *fd;
    int nb_read;
    double val;

    fprintf(stdout, "reading 'simple' vector b from '%s' ... ", filename);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }

    nb_read = fscanf(fd, "%ld", &n);
    if (nb_read != 1) {
        fprintf(stderr, "Error reading 'n' in '%s'\n", filename);
        exit(1);
    }
    fprintf(stderr, "n=%ld\n", n);

    b = new_vector(n);

    /* load vector data from file */
    j = 0;
    while (!feof(fd)) {
        nb_read = fscanf(fd, "%lf\n", &val);
        if (nb_read != 1) {
            fprintf(stderr, "Error reading mat[%ld] in '%s'\n",
                    j, filename);
            exit(1);
        }
        b->mat[j] = val;
        j++;
    }

    if (j != n) {
        fprintf(stderr,
                "read_simple_vector: truncated file ! read only %ld/%ld items.\n",
                j, n);
        exit(1);
    }

    fclose(fd);
    fprintf(stdout, "%ld lines\n", j);
    fflush(stdout);
    return (b);
}

/** \brief read a file as vector
 * 
 * the file is formated as follow :
 *
 * nb_total_of_item_in_vector
 * index1 value1
 * index2 value2
 * ...
 */
struct vector_t *read_vector(char *filename)
{
    struct vector_t *b;
    long int n, j;
    FILE *fd;
    int nb_read;
    long int rayid;
    double val;

    fprintf(stdout, "reading vector b from '%s' ... ", filename);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }

    nb_read = fscanf(fd, "%ld", &n);
    if (nb_read != 1) {
        fprintf(stderr, "Error reading 'n' in '%s'\n", filename);
        exit(1);
    }
    fprintf(stdout, "size is %ld ... ", n);
    fflush(stdout);

    b = new_vector(n);

    /* load vector data from file */
    j = 0;
    while (!feof(fd)) {
        nb_read = fscanf(fd, "%ld %lf\n", &rayid, &val);

        if (nb_read != 2) {
            fprintf(stderr, "Error reading mat[%ld] in '%s'\n",
                    j, filename);
            exit(1);
        }

        if (fabs(b->mat[rayid]) > 1.0e-6) {
            fprintf(stdout,
                    "read_vector: duplicate value (%ld) old=%f/new=%f\n",
                    rayid, b->mat[rayid], val);
        }
        b->mat[rayid] = val;
        j++;
    }
    fclose(fd);
    fprintf(stdout, "%ld lines\n", j);
    fflush(stdout);
    return (b);
}

void write_vector(struct vector_t *b, char *filename)
{
    long int i;
    FILE *fd;

    fprintf(stdout, "writing vector to '%s' ... ", filename);
    if (!(fd = fopen(filename, "w"))) {
        perror(filename);
        exit(1);
    }

    fprintf(fd, "%ld\n", b->length);
    for (i = 0; i < b->length; i++) {
        fprintf(fd, "%ld %f\n", i, b->mat[i]);
    }

    fclose(fd);
    fprintf(stdout, "%ld items\n", b->length);
}

/** \brief read a portion of a vector 
 *
 * very specific, used by ray2mesh to re-number the residuals.
 */
struct vector_t *read_subvector(char *filename, long int *first,
                                long int *last)
{
    struct vector_t *b;
    long int n, j;
    FILE *fd;
    int nb_read;
    long int rayid, min_rayid, max_rayid;
    double val;

    fprintf(stdout, "reading sub-vector b from '%s' ... ", filename);
    if (!(fd = fopen(filename, "r"))) {
        perror(filename);
        exit(1);
    }

    nb_read = fscanf(fd, "%ld", &n);
    if (nb_read != 1) {
        fprintf(stderr, "Error reading 'n' in '%s'\n", filename);
        exit(1);
    }

    b = new_vector(n);
    min_rayid = n;
    max_rayid = 0;

    /* load vector data from file */
    j = 0;
    while (1) {
        nb_read = fscanf(fd, "%ld %lf", &rayid, &val);

        if (feof(fd)) {
            break;
        }

        if (nb_read != 2) {
            fprintf(stderr, "Error reading mat[%ld] in '%s'\n",
                    j, filename);
            exit(1);
        }

        if (max_rayid < rayid) {
            max_rayid = rayid;
        }
        if (min_rayid > rayid) {
            min_rayid = rayid;
        }

        b->mat[rayid] = val;
        j++;
    }
    fclose(fd);

    *first = min_rayid;
    *last = max_rayid;

    b->length = max_rayid + 1;
    b->mat = (double *) realloc(b->mat, sizeof(double) * b->length);

    fprintf(stdout, "%ld lines\n", j);
    fflush(stdout);
    return (b);
}

struct vector_t *vector_resize(struct vector_t *v, long int new_length)
{
    long int old_length, i;

    if (!v) {
        return (NULL);
    }
    old_length = v->length;
    v->mat = (double *) realloc(v->mat, sizeof(double) * new_length);
    assert(v->mat);
    v->length = new_length;
    for (i = old_length; i < new_length; i++) {
        v->mat[i] = 0;
    }

    fprintf(stdout, "vector_resize: resize (%p) from %ld to %ld\n",
            v, old_length, new_length);
    return (v);
}

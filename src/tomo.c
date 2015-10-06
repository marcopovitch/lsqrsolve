#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <signal.h>
#include <libgen.h>

#include <mesh/mesh.h>
#include <mesh/layer.h>
#include <mesh/mesh2xml.h>
#include <mesh/import.h>
#include <mesh/metacell.h>

#include <ray/vitesse.h>
#include <sparse/sparse.h>

#ifdef __APPLE__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#include "lsqr.h"
#include "lsqr_wrapper.h"
#include "export.h"
#include "parse_options.h"
#include "irm.h"
#include "relative_tt.h"
#include "catch_sig.h"
#include "extern.h"
#include "regul.h"

/* declaration */
void check_files_access(int format, struct mesh_data_t *md, char *xmlfile);
void check_write_access(char *filename);
void sparse_compress_column(struct mesh_t *mesh, struct sparse_matrix_t *m,
                            long int **v);
struct vector_t *uncompress_column(struct vector_t *sol, long int *cmp2fat,
                                   long int nfat_sol);
void lsqrsolve_info();

/* tell if CTRL-C was pressed */
int please_stop_lsqr = 0;

/* dump intermediate solution */
int please_dump_lsqr = 0;

/********/
/* MAIN */
/********/
int main(int argc, char *argv[])
{
    struct sparse_matrix_t *sparseA = NULL;
    struct vector_t *b = NULL;
    struct vector_t *x;
    struct mesh_t *mesh;
    char *xml_output;

    long int *compress2fat = NULL;
    struct vector_t *solution;
    struct vector_t *std_error_sol;
    long int fat_sol_nb_col;

    lsqr_input *input;
    lsqr_output *output;
    lsqr_work *work;            /* zone temoraire de travail */
    lsqr_func *func;            /* func->mat_vec_prod -> APROD */

    /* cmd line arg */
    char *mesh_filename = NULL;
    char *importfilename = NULL;
    char *output_filename = NULL;
    char *sol_error_filename = NULL;
    char *log_filename = NULL;
    char *output_type = NULL;
    int max_iter;
    double damping, grad_damping;
    int use_ach = 0;            /* ACH : tele-seismic inversion tomography */
    int check_sparse = 0;       /* check sparse matrix disable by default  */

    /* velocity model */
    char *vmodel = NULL;
    struct velocity_model_t *vm = NULL;

    struct mesh_t **imported_mesh = NULL;
    char **xmlfilelist = NULL;
    int nb_xmlfile = 0;

    int i, j;

    int nb_irm = 0;
    struct irm_t **irm = NULL;
    int *nb_metacell = NULL;

    FILE *logfd;

    /*************************************************************/
    parse_command_line(argc, argv,
                       &mesh_filename,
                       &vmodel,
                       &importfilename,
                       &log_filename,
                       &output_filename, &output_type,
                       &max_iter,
                       &damping, &grad_damping, &use_ach, &check_sparse);

    if (use_ach) {
        fprintf(stderr, "Using ACH tomographic inversion\n");
    } else {
        fprintf(stderr, "Using STANDARD tomographic inversion\n");
    }

    /* load the velocity model */
    if (vmodel) {
        char *myfile;
        vm = load_velocity_model(vmodel);
        if (!vm) {
            fprintf(stderr, "Can not initialize velocity model '%s'\n",
                    vmodel);
            exit(1);
        }
        myfile = strdup(vmodel);
        fprintf(stderr, "Velocity model '%s' loaded\n", basename(myfile));
        free(myfile);
    } else {
        vm = NULL;
    }

    /* Open log file */
    if (!log_filename) {
        logfd = stdout;
    } else {
        if (!(logfd = fopen(log_filename, "w"))) {
            perror(log_filename);
            exit(1);
        }
    }

    /*check_write_access (output_filename); */

    /**************************************/
    /* test if we can open file to import */
    /**************************************/
    if (importfilename) {
        xmlfilelist = parse_separated_list(importfilename, ",");
        nb_xmlfile = 0;
        while (xmlfilelist[nb_xmlfile]) {
            if (access(xmlfilelist[nb_xmlfile], R_OK) == -1) {
                perror(xmlfilelist[nb_xmlfile]);
                exit(1);
            }
            nb_xmlfile++;
        }
    } else {
        fprintf(stderr, "No file to import ... exiting\n");
        exit(0);
    }

    /****************************/
    /* main mesh initialization */
    /****************************/
    mesh = mesh_init_from_file(mesh_filename);
    if (!mesh) {
        fprintf(stderr, "Error decoding %s.\n", mesh_filename);
        exit(1);
    }
    fprintf(stderr, "read %s ok\n", mesh_filename);

    /*****************************************/
    /* check and initialize slice xml files  */
    /*****************************************/
    if (nb_xmlfile) {
        int nb_sparse = 0;
        int nb_res = 0;
        int f;

        imported_mesh =
            (struct mesh_t **) malloc(sizeof(struct mesh_t *) *
                                      nb_xmlfile);
        assert(imported_mesh);

        for (i = 0; i < nb_xmlfile; i++) {
            imported_mesh[i] = mesh_init_from_file(xmlfilelist[i]);
            if (!imported_mesh[i]) {
                fprintf(stderr, "Error decoding %s.\n", mesh_filename);
                exit(1);
            }
            for (f = 0; f < NB_MESH_FILE_FORMAT; f++) {
                /* mandatory field : res, sparse, and irm if provided */
                if (f == RES || f == SPARSE || f == IRM) {
                    check_files_access(f, imported_mesh[i]->data[f],
                                       xmlfilelist[i]);
                }
            }
            if (imported_mesh[i]->data[SPARSE]) {
                nb_sparse += imported_mesh[i]->data[SPARSE]->ndatafile;
            }

            if (imported_mesh[i]->data[RES]) {
                nb_res += imported_mesh[i]->data[RES]->ndatafile;
            }

            if (imported_mesh[i]->data[IRM]) {
                nb_irm += imported_mesh[i]->data[IRM]->ndatafile;
            }

        }

        if (!nb_sparse || !nb_res) {
            fprintf(stderr, "Error no sparse or res file available !\n");
            exit(0);
        }
    }

    /*********************************************/
    /* read and import the sparse(s) matrix(ces) */
    /*********************************************/
    for (i = 0; i < nb_xmlfile; i++) {
        if (!imported_mesh[i]->data[SPARSE]) {
            continue;
        }

        for (j = 0; j < imported_mesh[i]->data[SPARSE]->ndatafile; j++) {
            sparseA = import_sparse_matrix(sparseA,
                                           imported_mesh[i]->data[SPARSE]->
                                           filename[j]);
        }
    }

    if (check_sparse) {
        if (check_sparse_matrix(sparseA)) {
            exit(1);
        }
    }

    /*sparse_compute_length(sparseA, "length1.txt"); */
    fat_sol_nb_col = sparseA->nb_col;
    show_sparse_stats(sparseA);

    /*********************************************/
    /* read and import the residual time vector  */
    /*********************************************/
    for (i = 0; i < nb_xmlfile; i++) {
        if (!imported_mesh[i]->data[RES]) {
            continue;
        }

        for (j = 0; j < imported_mesh[i]->data[RES]->ndatafile; j++) {
            b = import_vector(b, imported_mesh[i]->data[RES]->filename[j]);
        }
    }

    /*************************************************/
    /* check compatibility between matrix and vector */
    /*************************************************/
    if (sparseA->nb_line != b->length) {
        fprintf(stderr,
                "Error, check your matrix/vector sizes (%ld/%ld)\n",
                sparseA->nb_line, b->length);
        exit(1);
    }

    /********************/
    /* show memory used */
    /********************/
#ifdef __APPLE__
    {
        struct mstats memusage;

        memusage = mstats();
        fprintf(stderr, "Memory used: %.2f MBytes\n",
                (float) (memusage.bytes_used) / (1024. * 1024));
    }
#else
    {
        struct mallinfo m_info;

        m_info = mallinfo();
        fprintf(stderr, "Memory used: %.2f MBytes\n",
                (float) (m_info.uordblks +
                         m_info.usmblks) / (1024. * 1024.));
    }
#endif

    /**************************************/
    /* relative traveltime mode           */
    /**************************************/
    if (use_ach) {
        int nb_evt_imported = 0;

        for (i = 0; i < nb_xmlfile; i++) {
            if (!imported_mesh[i]->data[EVT]) {
                continue;
            }

            for (j = 0; j < imported_mesh[i]->data[EVT]->ndatafile; j++) {
                relative_tt(sparseA, b,
                            imported_mesh[i]->data[EVT]->filename[j]);
                nb_evt_imported++;
            }
        }

        if (!nb_evt_imported) {
            fprintf(stderr,
                    "Error in ACH mode, can not import any .evt file !\n");
            exit(1);
        }
    }

    /************************************************/
    /* read the irregular mesh definition if needed */
    /* one by layer                                 */
    /************************************************/
    if (nb_irm) {
        int cpt = 0;
        struct mesh_offset_t **offset;
        int l;

        irm = (struct irm_t **) malloc(nb_irm * sizeof(struct irm_t *));
        assert(irm);
        nb_metacell = (int *) calloc(nb_irm, sizeof(int));
        assert(nb_metacell);

        make_mesh(mesh);

        for (i = 0; i < nb_xmlfile; i++) {
            if (!imported_mesh[i]->data[IRM]) {
                continue;
            }

            /* offset between meshes */
            offset = compute_mesh_offset(mesh, imported_mesh[i]);
            for (l = 0; l < mesh->nlayers; l++) {
                if (!offset[l])
                    continue;
                fprintf(stderr,
                        "\t%s, [%s] offset[layer=%d] : lat=%d lon=%d z=%d\n",
                        xmlfilelist[i], MESH_FILE_FORMAT[IRM], l,
                        offset[l]->lat, offset[l]->lon, offset[l]->z);
            }

            for (j = 0; j < imported_mesh[i]->data[IRM]->ndatafile; j++) {
                /* FIXME: read only once the irm file */
                irm[cpt] =
                    read_irm(imported_mesh[i]->data[IRM]->filename[j],
                             &(nb_metacell[cpt]));
                import2mesh_irm_file(mesh,
                                     imported_mesh[i]->data[IRM]->
                                     filename[j], offset);
                cpt++;
            }

            for (l = 0; l < mesh->nlayers; l++) {
                if (offset[l])
                    free(offset[l]);
            }
            free(offset);
        }
        metacell_find_neighbourhood(mesh);
    }

    /*sparse_compute_length(sparseA, "length1.txt"); */
    fat_sol_nb_col = sparseA->nb_col;
    show_sparse_stats(sparseA);

    /***********************/
    /* remove empty column */
    /***********************/
    fprintf(stderr, "starting compression ...\n");
    sparse_compress_column(mesh, sparseA, &compress2fat);
    if (check_sparse) {
        if (check_sparse_matrix(sparseA)) {
            exit(1);
        }
    }
    show_sparse_stats(sparseA);

    /***************************************/
    /* add gradient damping regularisation */
    /***************************************/
    if (fabs(grad_damping) > 1.e-6) {
        int nb_faces = 6;       /* 1 cell may have 6 neighbours */
        long int nb_lines = 0;
        char *regul_name;

        fprintf(stdout, "using gradient damping : %f\n", grad_damping);

        /* tmp file name */
        regul_name = tempnam("/tmp", "regul");
        if (!regul_name) {
            perror("lsqrsolve: ");
            exit(1);
        }

        if (nb_irm) {
            create_regul_DtD_irm(sparseA, compress2fat,
                                 mesh, regul_name, nb_faces,
                                 grad_damping, &nb_lines);
        } else {
            create_regul_DtD(sparseA, compress2fat,
                             mesh, regul_name, nb_faces,
                             grad_damping, &nb_lines);
        }

        sparse_matrix_resize(sparseA,
                             sparseA->nb_line + sparseA->nb_col,
                             sparseA->nb_col);

        sparseA = import_sparse_matrix(sparseA, regul_name);
        if (check_sparse) {
            if (check_sparse_matrix(sparseA)) {
                exit(1);
            }
        }
        vector_resize(b, sparseA->nb_line);
        unlink(regul_name);

        show_sparse_stats(sparseA);
    }

    /*********************************/
    /* the real mesh is no more used */
    /* keep only the light mesh      */
    /*********************************/
    fprintf(stdout,
            "Time to free the real mesh and keep only the light structure\n");
    free_mesh(mesh);
    mesh = mesh_init_from_file(mesh_filename);
    if (!mesh) {
        fprintf(stderr, "Error decoding %s.\n", mesh_filename);
        exit(1);
    }
    fprintf(stderr, "read %s ok\n", mesh_filename);

    /********************************/
    /* init vector solution to zero */
    /********************************/
    x = new_vector(sparseA->nb_col);

    /*************************************************************/
    /* solve A.x = B                                             */
    /* A = ray length in the cells                               */
    /* B = residual travel time observed - computed              */
    /* x solution to satisfy the lsqr problem                    */
    /*************************************************************/

    /* LSQR alloc */
    alloc_lsqr_mem(&input, &output, &work, &func, sparseA->nb_line,
                   sparseA->nb_col);

    fprintf(stderr, "alloc_lsqr_mem : ok\n");

    /* defines the routine Mat.Vect to use */
    func->mat_vec_prod = sparseMATRIXxVECTOR;

    /* Set the input parameters for LSQR */
    input->num_rows = sparseA->nb_line;
    input->num_cols = sparseA->nb_col;
    input->rel_mat_err = 1.0e-3;        /* in km */
    input->rel_rhs_err = 1.0e-2;        /* in seconde */
    /*input->rel_mat_err = 0.;
       input->rel_rhs_err = 0.; */
    input->cond_lim = .0;
    input->lsqr_fp_out = logfd;
    /* input->rhs_vec = (dvec *) b; */
    dvec_copy((dvec *) b, input->rhs_vec);
    input->sol_vec = (dvec *) x;        /* initial guess */
    input->damp_val = damping;
    if (max_iter == -1) {
        input->max_iter = 4 * (sparseA->nb_col);
    } else {
        input->max_iter = max_iter;
    }

    /* catch Ctrl-C signal */
    signal(SIGINT, emergency_halt);

    /******************************/
    /* resolution du systeme Ax=b */

    /******************************/
    lsqr(input, output, work, func, sparseA);
    fprintf(stderr, "*** lsqr ended (%ld iter) : %s\n",
            output->num_iters, lsqr_msg[output->term_flag]);
    if (output->term_flag == 0) {       /* solution x=x0 */
        exit(0);
    }

    /* uncompress the solution  */
    solution = uncompress_column((struct vector_t *) output->sol_vec,
                                 compress2fat, fat_sol_nb_col);

    /* uncompress the standard error on solution  */
    std_error_sol =
        uncompress_column((struct vector_t *) output->std_err_vec,
                          compress2fat, fat_sol_nb_col);

    /* if irm file was provided, set the right value to each cell 
     * from a given metacell 
     */
    if (irm) {
        irm_update(solution, irm, nb_metacell, nb_irm, mesh);
        free_irm(irm, nb_irm);
        free(nb_metacell);
    }

    /* write solution */
    if (strchr(output_type, 'm')) {
        export2matlab(solution, output_filename, mesh, vm,
                      output->num_iters,
                      input->damp_val, grad_damping, use_ach);
    }

    if (strchr(output_type, 's')) {
        export2sco(solution, output_filename, mesh, vm,
                   output->num_iters,
                   input->damp_val, grad_damping, use_ach);
    }

    if (strchr(output_type, 'g')) {
        /* solution */
        export2gmt(solution, output_filename, mesh, vm,
                   output->num_iters,
                   input->damp_val, grad_damping, use_ach);

        /* error on solution */
        sol_error_filename = (char *)
            malloc(sizeof(char) *
                   (strlen(output_filename) + strlen(".err") + 1));
        sprintf(sol_error_filename, "%s.err", output_filename);
        export2gmt(std_error_sol, sol_error_filename, mesh, vm,
                   output->num_iters,
                   input->damp_val, grad_damping, use_ach);
        free(sol_error_filename);

    }

    /* save the xml enrichied with sections */
    xml_output = (char *)
        malloc((strlen(output_filename) + strlen(".xml") +
                1) * sizeof(char));
    assert(xml_output);
    sprintf(xml_output, "%s.xml", output_filename);
    mesh2xml(mesh, xml_output);
    free(xml_output);

    /******************************************************/
    /* variance reduction, ie how the model fits the data */
    /* X = the final solution                             */
    /*                                                    */
    /*                 ||b-AX||²                          */
    /*         VR= 1 - --------                           */
    /*                  ||b||²                            */
    /*                                                    */
    /******************************************************/
    {
        double norm_b;
        double norm_b_AX;
        double VR;              /* variance reduction */

        struct vector_t *rhs;   /* right hand side */
        rhs = new_vector(sparseA->nb_line);

        /* use copy */
        dvec_copy((dvec *) b, (dvec *) rhs);

        norm_b = dvec_norm2((dvec *) rhs);

        /* does rhs = rhs + sparseA . output->sol_vec */
        /* here  rhs is overwritten */
        dvec_scale((-1.0), (dvec *) rhs);
        sparseMATRIXxVECTOR(0, output->sol_vec, (dvec *) rhs, sparseA);
        dvec_scale((-1.0), (dvec *) rhs);

        norm_b_AX = dvec_norm2((dvec *) rhs);

        VR = 1 - (norm_b_AX * norm_b_AX) / (norm_b * norm_b);
        fprintf(stdout, "Variance reduction = %.2f%%\n", VR * 100);
        free_vector(rhs);
    }

    /********/
    /* free */
    /********/
    if (vm) {
        free_velocity_model(vm);
    }
    free_mesh(mesh);
    free_sparse_matrix(sparseA);
    free_lsqr_mem(input, output, work, func);

    free_vector(solution);
    free_vector(std_error_sol);
    free(compress2fat);

    for (i = 0; i < nb_xmlfile; i++) {
        free(xmlfilelist[i]);
        free_mesh(imported_mesh[i]);
    }
    free(xmlfilelist);
    free(imported_mesh);
    return (0);
}





/** \brief remove all empty column witch does not carry head cell 
 *  from the sparse matrix 
 *
 * all empties cells (no ray/ no metacell) are virtually removed
 * from the inversion. Those cells will not
 * have any associated velocity in the output files.
 * This flag will also modified the regularization
 * behaviour, such as gradient smoothing is only done
 * with non empty cells.
 * 
 * return a new smaller sparse matrix and a vector to keep the real
 * rank of each column.
 */
void sparse_compress_column(struct mesh_t *mesh, struct sparse_matrix_t *m,
                            long int **v)
{
    long int j;
    long int nb;
    long int *vect;
    struct sparse_item_t *cur_item;
    struct sparse_item_t **new_col;
    struct sparse_item_t **new_last_col;
    int x, y, z;
    struct cell_t *hcell;

    if (mesh->nb_total_metacell) {
        fprintf(stdout, "%s: looking for null column without metacell\n",
                __FUNCTION__);
        assert(mesh->cell);

    } else {
        fprintf(stdout, "%s: looking for null column\n", __FUNCTION__);
    }

    /* keep trace of non empty column */
    vect = (long int *) calloc(m->nb_col, sizeof(long int));
    assert(vect);

    /* new compressed columns */
    new_col = (struct sparse_item_t **)
        calloc(m->nb_col, sizeof(struct sparse_item_t *));
    assert(new_col);

    /* copy the "last_col" info */
    new_last_col = (struct sparse_item_t **)
        calloc(m->nb_col, sizeof(struct sparse_item_t *));

    /* how many non-empty column */
    nb = 0;
    for (j = 0; j < m->nb_col; j++) {

        if (!m->col[j]) {
            if (!mesh->nb_total_metacell) {
                continue;
            }
            /* check this cell is not a head cell */
            unlinearize_cell_id(j, &x, &y, &z, mesh);
            hcell = move_in_mesh(mesh, NULL, x, y, z);
            if (!hcell->nb_meta_scores) {
                continue;
            }
        }

        /* non empty cell : keep it */
        vect[nb] = j;
        new_last_col[nb] = m->last_col[j];
        new_col[nb] = m->col[j];

        /* column renumbering of sparse items */
        cur_item = m->col[j];
        while (cur_item) {
            cur_item->col_index = nb;
            cur_item = cur_item->next_in_col;
        }

        nb++;
    }

    /* resize */
    vect = (long int *) realloc(vect, sizeof(long int) * nb);
    assert(vect);
    *v = vect;

    new_col = (struct sparse_item_t **)
        realloc(new_col, sizeof(struct sparse_item_t *) * nb);
    assert(new_col);

    new_last_col = (struct sparse_item_t **)
        realloc(new_last_col, sizeof(struct sparse_item_t *) * nb);
    assert(new_last_col);

    /* free */
    free(m->col);
    m->col = new_col;
    m->nb_col = nb;

    free(m->last_col);
    m->last_col = new_last_col;
}

struct vector_t *uncompress_column(struct vector_t *sol, long int *cmp2fat,
                                   long int nfat_sol)
{
    long int i;
    struct vector_t *fat_sol;

    fat_sol = new_vector(nfat_sol);
    assert(fat_sol);

    /* set to "novalue" the cells wich have not been processed by lsqr */
    for (i = 0; i < nfat_sol; i++) {
        fat_sol->mat[i] = NOVALUE;
    }

    for (i = 0; i < sol->length; i++) {
        fat_sol->mat[cmp2fat[i]] = sol->mat[i];
    }

    return (fat_sol);
}

void check_files_access(int format, struct mesh_data_t *md, char *xmlfile)
{
    int i;
    char *format_str;

    if (!md) {
        format_str = mesh_get_section_name(format);
        fprintf(stderr,
                "warning: %s does not have any associated %s file(s)\n",
                xmlfile, format_str);
        free(format_str);
        return;
    }

    for (i = 0; i < md->ndatafile; i++) {
        if (access(md->filename[i], R_OK) == -1) {
            perror(md->filename[i]);
            exit(1);
        }
    }
}

void check_write_access(char *filename)
{
    if (access(filename, W_OK) < 0) {
        fprintf(stderr,
                "%s : cannot access file %s (write mode). Exiting.\n",
                "solve", filename);
        exit(1);
    }
}

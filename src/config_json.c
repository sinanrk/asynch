#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <memory.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(_MSC_VER)
#include <process.h>
#endif

#include <structs.h>
#include <globals.h>
#include <config_json.h>
#include <cson.h>
#include <db.h>
#include <outputs.h>
#include <io.h>
#include <models/definitions.h>


#if defined(_MSC_VER)
#   define timegm _mkgmtime
#endif


//Print to stderr only for process of rank 0
static int print_err(const char* format, ...)
{
    int res = 0;
    if (my_rank == 0)
    {
        va_list args;
        va_start(args, format);
        res = vfprintf(stderr, format, args);
        va_end(args);
    }

    return res;
}

// Getting file extension
static const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if (!dot || dot == filename) return "";
    return dot + 1;
}

#define CSON_GET_VALUE(obj, key, v) \
v = cson_object_get((obj), (key)); \
if (v == NULL) \
{ \
    print_err("Error: Did not get a value from config file for %s.\n", (key)); \
    return NULL; \
}

#define CSON_GET_SUB_VALUE(obj, key, v) \
v = cson_object_get_sub((obj), (key), '.'); \
if (v == NULL) \
{ \
    print_err("Error: Did not get a value from config file for %s.\n", (key)); \
    return NULL; \
}

#define CSON_CHECK_VALUE(key, type, v) \
if (!cson_value_is_ ## type ## (v)) \
{ \
    print_err("Error: Invalid " key " type (" #type " expected).\n"); \
    return NULL; \
}

GlobalVars* Read_Config_JSON(
    const char * const filename,
    ErrorData *errors,
    Forcing *forcings,
    ConnData *db_connections,
    char *rkdfilename,
    AsynchModel const *model,
    void *external)
{
    GlobalVars* globals = malloc(sizeof(GlobalVars));
    memset(globals, 0, sizeof(GlobalVars));

    FILE *in = fopen(filename, "r");
    if (in == NULL)
    {
        printf("Error: Global file %s was not found.\n", filename);
        return NULL;
    }

    cson_value * root = NULL;
    int rc = cson_parse_FILE(&root, in, NULL, NULL); 

    // The NULL arguments hold optional information for/about
    // the parse results. These can be used to set certain
    // parsing options and get more detailed error information
    // if parsing fails.

    if (0 != rc)
    {
        printf("Error: Failed to parse config file (%s).\n", cson_rc_string(rc));
        return NULL;
    }

    if (!cson_value_is_object(root))
    {
        print_err("Error: Invalid config file (root object expected).\n");
        return NULL;
    }

    cson_object * obj = cson_value_get_object(root);

    //Grab the model uid
    cson_value * v;
    CSON_GET_VALUE(obj, "model", v);
    CSON_CHECK_VALUE("model", integer, v);
    globals->model_uid = (unsigned short) cson_value_get_integer(v);

    //Grab the begin and end time
    struct tm begin_tm;
    memset(&begin_tm, 0, sizeof(struct tm));

    CSON_GET_VALUE(obj, "begin", v);
    if (cson_value_is_integer(v))
        globals->begin_time = cson_value_get_integer(v);
    else
    {
        int valsread = sscanf(cson_value_get_cstr(v), "%d-%d-%d %d:%d", &begin_tm.tm_year, &begin_tm.tm_mon, &begin_tm.tm_mday, &begin_tm.tm_hour, &begin_tm.tm_min);
        if (valsread == 5)
        {
            begin_tm.tm_year = begin_tm.tm_year - 1900;
            begin_tm.tm_mon = begin_tm.tm_mon - 1;
            globals->begin_time = timegm(&begin_tm);
        }
        else
        {
            print_err("Error: Invalid begin datetime format (expected YYYY-MM-DD HH:MM || unix_time, got %s).\n", cson_value_get_cstr(v));
            return NULL;
        }
    }

    struct tm end_tm;
    memset(&end_tm, 0, sizeof(struct tm));

    CSON_GET_VALUE(obj, "end", v);
    if (cson_value_is_integer(v))
        globals->end_time = cson_value_get_integer(v);
    else
    {
        int valsread = sscanf(cson_value_get_cstr(v), "%d-%d-%d %d:%d", &end_tm.tm_year, &end_tm.tm_mon, &end_tm.tm_mday, &end_tm.tm_hour, &end_tm.tm_min);
        if (valsread == 5)
        {
            end_tm.tm_year = end_tm.tm_year - 1900;
            end_tm.tm_mon = end_tm.tm_mon - 1;
            globals->end_time = timegm(&end_tm);
        }
        else
        { 
            print_err("Error: Invalid end datetime format (expected YYYY-MM-DD HH:MM || unix_time, got %s).\n", cson_value_get_cstr(v));
            return NULL;
        }
    }

    globals->maxtime = (double)(globals->end_time - globals->begin_time) / 60.0;
    if (globals->maxtime <= 0.0)
    {
        print_err("Error: Simulation period invalid (begin >= end)\n");
        return NULL;
    }

    //Grab the outputs
    cson_value *outputs_obj = cson_object_get(obj, "outputs");
    if (!cson_value_is_object(outputs_obj))
    {
        print_err("Error: Invalid outputs (object expected).\n");
        return NULL;
    }

    {
        cson_object *obj = cson_value_get_object(outputs_obj);

        // Default to false
        globals->with_output_postfix_param = false;

        v = cson_object_get(obj, "postfix_with_global_params");
        if (v && cson_value_is_bool(v))
            globals->with_output_postfix_param = cson_value_get_bool(v);

        cson_value *variables = cson_object_get(obj, "variables");
        if (cson_value_is_array(variables))
        {
            cson_array *ar = cson_value_get_array(variables);
            assert(ar != NULL);

            globals->num_outputs = cson_array_length_get(ar);
            if (globals->num_outputs > 0)
            {
                globals->outputs = malloc(globals->num_outputs * sizeof(Output));

                for (unsigned int i = 0; i < globals->num_outputs; ++i)
                {
                    v = cson_array_get(ar, i);
                    globals->outputs[i].name = strdup(cson_value_get_cstr(v));
                }
            }
        }
        else
        {
            print_err("Error: Invalid outputs.variables (expected array)");
            return NULL;
        }

        if (cson_object_get(obj, "timeseries"))
        {
            //Grab where to write the timeseries
            CSON_GET_SUB_VALUE(obj, "timeseries.filename", v);
            globals->hydros_loc_filename = strdup(cson_value_get_cstr(v));

            const char *ext = get_filename_ext(globals->hydros_loc_filename);
            if (strcmp(ext, "dat") == 0)
                globals->hydros_loc_flag = 1;
            else if (strcmp(ext, "csv") == 0)
                globals->hydros_loc_flag = 2;
            else if (strcmp(ext, "dbc") == 0)
            {
                globals->hydros_loc_flag = 3;

                // TODO
            }                
            else if (strcmp(ext, "rad") == 0)
                globals->hydros_loc_flag = 4;
            else if (strcmp(ext, "h5") == 0)
                globals->hydros_loc_flag = 5;
            else
            {
                print_err("Error: Invalid timeseries.filename extension (expected .dat, .csv, .dbc or .h5)");
                return NULL;
            }

            v = cson_object_get_sub(obj, "timeseries.locations", '.');
            if (v)
            {
                if (cson_value_is_string(v))
                {
                    globals->hydrosave_filename = strdup(cson_value_get_cstr(v));

                    const char *ext = get_filename_ext(globals->hydrosave_filename);
                    if (strcmp(ext, "sav") == 0)
                        globals->hydrosave_flag = 1;
                    else if (strcmp(ext, "dbc") == 0)
                    {
                        globals->hydrosave_flag = 3;

                        // TODO
                    }
                    else
                    {
                        print_err("Error: Invalid timeseries.locations file extension (expected .sav or .dbc)");
                        return NULL;
                    }
                }
                else
                {
                    print_err("Error: Invalid timeseries.locations location (expected string)");
                    return NULL;
                }
            }
            else
                // Default to all links
                globals->hydrosave_flag = 3;

            CSON_GET_SUB_VALUE(obj, "timeseries.interval", v);
            globals->print_time = cson_value_get_double(v);
        }

        if (cson_object_get(obj, "peaks"))
        {
            //Grab where to write the timeseries
            CSON_GET_SUB_VALUE(obj, "peaks.filename", v);
            globals->peaks_loc_filename = strdup(cson_value_get_cstr(v));

            const char *ext = get_filename_ext(globals->peaks_loc_filename);
            if (strcmp(ext, "pea") == 0)
                globals->peaks_loc_flag = 1;
            else if (strcmp(ext, "dbc") == 0)
            {
                globals->peaks_loc_flag = 2;

                // TODO .dbc
            }
            else
            {
                print_err("Error: Invalid peaks file extension (expected .pea or .dbc)");
                return NULL;
            }

            v = cson_object_get_sub(obj, "peaks.locations", '.');
            if (v)
            {
                if (cson_value_is_string(v))
                {
                    globals->peaksave_filename = strdup(cson_value_get_cstr(v));

                    const char *ext = get_filename_ext(globals->peaksave_filename);
                    if (strcmp(ext, "sav") == 0)
                        globals->peaksave_flag = 1;
                    else if (strcmp(ext, "dbc") == 0)
                    {
                        globals->peaksave_flag = 3;

                        // TODO .dbc
                    }
                    else
                    {
                        print_err("Error: Invalid peaks location file extension (expected .sav or .dbc)");
                        return NULL;
                    }
                }
                else
                {
                    print_err("Error: Invalid peaks location (expected string)");
                    return NULL;
                }
            }
            else
                // Default to all links
                globals->hydrosave_flag = 3;

            v = cson_object_get_sub(obj, "peaks.function", '.');
            if (v)
            {
                globals->peakflow_function_name = strdup(cson_value_get_cstr(v));
                SetPeakflowOutputFunctions(globals->peakflow_function_name, &globals->peakflow_output);
            }
            else
                // Default to Classic
                globals->peakflow_output = &OutputPeakflow_Classic_Format;
        }

    }

    // Global parameters

    //Grab the global parameters
    //ReadLineFromTextFile(globalfile, line_buffer, line_buffer_len);
    //valsread = sscanf(line_buffer, "%u%n", &globals->num_global_params, &total);
    //if (ReadLineError(valsread, 1, "number of global parameters"))	return NULL;
    //globals->global_params = malloc(globals->num_global_params * sizeof(double));
    ////if (my_rank == 0 && globals->num_global_params != model->num_global_params)
    ////{
    ////    printf("[%i]: Error: Got %u global params in the .gbl file. Expected %u for model %u.\n", my_rank, globals->num_global_params, model->num_global_params, globals->model_uid);
    ////    MPI_Abort(MPI_COMM_WORLD, 1);
    ////}
    //for (i = 0; i < globals->num_global_params; i++)
    //{
    //    valsread = sscanf(&(line_buffer[total]), "%lf%n", &globals->global_params[i], &written);
    //    if (ReadLineError(valsread, 1, "a global parameter"))	return NULL;
    //    total += written;
    //}

    cson_value *global_params_obj = cson_object_get(obj, "global_params");
    if (cson_value_is_array(global_params_obj))
    {
        cson_array *ar = cson_value_get_array(global_params_obj);
        assert(ar != NULL);

        globals->num_global_params = cson_array_length_get(ar);
        if (globals->num_global_params > 0)
        {
            globals->global_params = malloc(globals->num_global_params * sizeof(double));

            for (unsigned int i = 0; i < globals->num_global_params; ++i)
            {
                v = cson_array_get(ar, i);
                globals->global_params[i] = cson_value_get_double(v);
            }
        }
    }

    //Set dim and other sizes
    if (model)
        model->set_param_sizes(globals, external);
    else
        SetParamSizes(globals, external);

    //Find the states needed for printing
    globals->num_states_for_printing = 0;
    globals->print_indices = calloc(globals->num_outputs, sizeof(unsigned int));
    for (unsigned int i = 0; i < globals->num_outputs; i++)
        SetDefaultOutputFunctions(globals->outputs[i].name, &globals->outputs[i], globals->print_indices, &globals->num_states_for_printing);

    globals->print_indices = realloc(globals->print_indices, globals->num_states_for_printing * sizeof(unsigned int));

    //Grab the stored steps limits
    v = cson_object_get_sub(obj, "buffers.num_step", '.');
    if (v)
    {
        CSON_CHECK_VALUE("buffers.num_step", integer, v);
        globals->iter_limit = (int) cson_value_get_integer(v);
    }        
    else 
        globals->iter_limit = 30;

    v = cson_object_get_sub(obj, "buffers.num_transfer", '.');
    if (v)
    {
        CSON_CHECK_VALUE("buffers.num_transfer", integer, v);
        globals->max_transfer_steps = (int) cson_value_get_integer(v);
    }        
    else
        globals->max_transfer_steps = 10;

    v = cson_object_get_sub(obj, "buffers.num_discont", '.');
    if (v)
    {
        CSON_CHECK_VALUE("buffers.num_discont", integer, v);
        globals->discont_size = (unsigned int) cson_value_get_integer(v);
    }        
    else
        globals->discont_size = 10;

    //Grab the topology data filename
    {
        CSON_GET_VALUE(obj, "topology", v);
        CSON_CHECK_VALUE("topology", string, v);
        globals->rvr_filename = strdup(cson_value_get_cstr(v));

        const char *ext = get_filename_ext(globals->rvr_filename);
        if (strcmp(ext, "rvr") == 0)
            globals->rvr_flag = 0;
        else if (strcmp(ext, "dbc") == 0)
        {
            globals->rvr_flag = 1;

            // TODO .dbc
        }
        else
        {
            print_err("Error: Invalid topology file extension (expected .rvr or .dbc)");
            return NULL;
        }
    }

    //Grab the local parameter data filename
    {
        CSON_GET_VALUE(obj, "local_params", v);
        CSON_CHECK_VALUE("local_params", string, v);
        globals->prm_filename = strdup(cson_value_get_cstr(v));

        const char *ext = get_filename_ext(globals->prm_filename);
        if (strcmp(ext, "prm") == 0)
            globals->prm_flag = 0;
        else if (strcmp(ext, "dbc") == 0)
        {
            globals->prm_flag = 1;

            // TODO .dbc
        }
        else
        {
            print_err("Error: Invalid local parameter file extension (expected .prm or .dbc)");
            return NULL;
        }
    }

    //Grab the initial state data filename
    {
        CSON_GET_VALUE(obj, "initial_state", v);
        CSON_CHECK_VALUE("initial_state", string, v);
        globals->init_filename = strdup(cson_value_get_cstr(v));

        const char *ext = get_filename_ext(globals->init_filename);
        if (strcmp(ext, "ini") == 0)
            globals->init_flag = 0;
        else if (strcmp(ext, "uini") == 0)
            globals->init_flag = 1;
        else if (strcmp(ext, "rec") == 0)
            globals->init_flag = 2;
        else if (strcmp(ext, "dbc") == 0)
        {
            globals->init_flag = 3;

            // TODO .dbc
        }
        else if (strcmp(ext, "h5") == 0)
            globals->init_flag = 4;
        else
        {
            print_err("Error: Invalid local parameter file extension (expected .prm or .dbc)");
            return NULL;
        }
    }

    //Grab forcings
    cson_value *forcings_obj = cson_object_get(obj, "forcings");
    if (!cson_value_is_array(forcings_obj))
    {
        print_err("Error: Invalid forcings (object expected).\n");
        return NULL;
    }

    {
        cson_array *ar = cson_value_get_array(forcings_obj);
        assert(ar != NULL);

        unsigned int num_forcings = cson_array_length_get(ar);        
        if (num_forcings < globals->num_forcings)
        {
            print_err("Error: Got %u forcings in the .gbl file (expected %u for model %u).\n", num_forcings, globals->num_forcings, globals->model_uid);
            return NULL;
        }
        if (num_forcings > globals->num_forcings && my_rank == 0)
        {
            print_err("Warning: Got %u forcings in the .gbl file. Expected %u for model %u.\n", my_rank, num_forcings, globals->num_forcings, globals->model_uid);
            globals->num_forcings = num_forcings;
        }
        
        for (unsigned int i = 0; i < num_forcings; ++i)
        {
            v = cson_array_get(ar, i);

            if (cson_value_is_null(v))
                forcings[i].flag = 0;
            if (cson_value_is_string(v))
            {
                forcings[i].filename = strdup(cson_value_get_cstr(v));

                const char *ext = get_filename_ext(forcings[i].filename);                
                if (strcmp(ext, "str") == 0)
                    globals->init_flag = 1;
                else if (strcmp(ext, "ustr") == 0)
                    globals->init_flag = 4;
                else if (strcmp(ext, "mon") == 0)
                {
                    globals->init_flag = 7;
                    forcings[i].raindb_start_time = globals->begin_time;
                    forcings[i].first_file = globals->begin_time;
                    forcings[i].last_file = globals->end_time;
                }
                else
                {
                    print_err("Error: Invalid forcing file extension (expected .str or .ustr)");
                    return NULL;
                }

            }
            else if (cson_value_is_object(v))
            {
                cson_object * forcing_obj = cson_value_get_object(v);

                CSON_GET_VALUE(forcing_obj, "filename", v);
                CSON_CHECK_VALUE("forcing.filename", string, v);

                forcings[i].filename = strdup(cson_value_get_cstr(v));

                const char *ext = get_filename_ext(forcings[i].filename);
                if (strcmp(ext, "bin") == 0)
                {
                    globals->init_flag = 3;
                }
                else if (strcmp(ext, "dbc") == 0)
                {
                    globals->init_flag = 3;

                    CSON_GET_VALUE(forcing_obj, "chunk_size", v);
                    CSON_CHECK_VALUE("forcing.chunk_size", integer, v);
                    forcings[i].increment = (unsigned int) cson_value_get_integer(v); 

                    CSON_GET_VALUE(forcing_obj, "time_step", v);
                    CSON_CHECK_VALUE("forcing.time_step", double, v);
                    forcings[i].file_time = cson_value_get_double(v);

                    forcings[i].raindb_start_time = globals->begin_time;
                    forcings[i].first_file = globals->begin_time;
                    forcings[i].last_file = globals->end_time;
                }
                else if (strcmp(ext, "gz") == 0)
                {
                    globals->init_flag = 6;

                    CSON_GET_VALUE(forcing_obj, "chunk_size", v);
                    CSON_CHECK_VALUE("forcing.chunk_size", integer, v);
                    forcings[i].increment = (unsigned int)cson_value_get_integer(v);

                    CSON_GET_VALUE(forcing_obj, "time_step", v);
                    CSON_CHECK_VALUE("forcing.time_step", double, v);
                    forcings[i].file_time = cson_value_get_double(v);
                }
                else if (strcmp(ext, "mon") == 0)
                {
                    globals->init_flag = 7;

                    forcings[i].raindb_start_time = globals->begin_time;
                    forcings[i].first_file = globals->begin_time;
                    forcings[i].last_file = globals->end_time;
                }
                else if (strcmp(ext, "cel") == 0)
                {
                    globals->init_flag = 8;

                    CSON_GET_VALUE(forcing_obj, "chunk_size", v);
                    CSON_CHECK_VALUE("forcing.chunk_size", integer, v);
                    forcings[i].increment = (unsigned int)cson_value_get_integer(v);

                    CSON_GET_VALUE(forcing_obj, "first_file", v);
                    CSON_CHECK_VALUE("forcing.first_file", double, v);
                    forcings[i].first_file = (unsigned int)cson_value_get_integer(v);

                    CSON_GET_VALUE(forcing_obj, "last_file", v);
                    CSON_CHECK_VALUE("forcing.last_file", double, v);
                    forcings[i].last_file = (unsigned int)cson_value_get_integer(v);
                }
                else
                {
                    print_err("Error: Invalid forcing file extension (expected .bin, .dbc, .gz, .mon or .cel)");
                    return NULL;
                }

            }
            else
            {
                print_err("Error: Invalid forcing (string or object expected).\n");
                return NULL;
            }


        }
    }

    // Clean up
    cson_value_free(root);

    return globals;
}

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#include <rksteppers.h>
#include <models/definitions.h>
#include <models/equations.h>
#include <models/check_consistency.h>
#include <models/output_constraints.h>
#include <models/check_state.h>

//Sets the various sizes and flags for the model. This method should set the following fields:
//dim:			The number of unknowns in the differential equations (or the number of ODEs at each link).
//diff_start:		The starting index of the differential unknowns.
//str_flag:		1 if reading rainfall data from a .str file, 0 else.
//binrain_flag:		1 if reading rainfall data from binary files, 0 else.
//uses_dam:		1 if dams are compatible with the model given by model_uid.
//params_size:		The total number of parameters (including precalculations) at links with no dam.
//dam_params_size:	The total number of parameters (including precalculations) at links with a dam.
//area_idx:		The entry in params where the upstream area is stored.
//disk_params:		The number of enries in param that are read from DEM data.
//Currently, this program assumes the same number of differential equations at each link.
//UnivVars* GlobalVars:	Contains the global variables for the system.
void SetParamSizes(GlobalVars* globals, void* external) {
	unsigned short int model_uid = globals->model_uid;
	unsigned int num_global_params;

    //Set dim and start of differential variables
    switch (model_uid)
    {
        //--------------------------------------------------------------------------------------------
    case 0:	num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 1;
        globals->num_forcings = 0;
        globals->min_error_tolerances = 1;
        break;
        //--------------------------------------------------------------------------------------------
//    case 1: num_global_params = 6;
//        globals->uses_dam = 0;
//        globals->num_params = 20;
//        globals->dam_params_size = 0;
//        globals->area_idx = 2;
//        globals->areah_idx = 1;
//        globals->num_disk_params = 12;
//        globals->convertarea_flag = 1;
//        globals->num_forcings = 1;
//        globals->min_error_tolerances = 1;	//This should probably be higher...
//        break;

        case 1: num_global_params = 6;
            globals->uses_dam = 0;
            globals->num_params = 20;
            globals->dam_params_size = 0;
            globals->area_idx = 2;
            globals->areah_idx = 1;
            globals->num_disk_params = 12;
            globals->convertarea_flag = 1;
            globals->num_forcings = 1;
            globals->min_error_tolerances = 1;	//This should probably be higher...
            break;
        //--------------------------------------------------------------------------------------------
    case 2: num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 3: num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 4: num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 5: num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 6: num_global_params = 5;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 14;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 15: 	num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 19:	num_global_params = 7;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
  
case 20:	num_global_params = 9;
        globals->uses_dam = 0;
        globals->num_params = 6;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 21:	num_global_params = 7;
        globals->uses_dam = 1;
        globals->num_params = 7;	//Need 1 extra for orifice_area
        globals->dam_params_size = 15;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 22:	num_global_params = 4;
        globals->uses_dam = 1;
        globals->num_params = 10;		//Need 1 extra for orifice_area
        globals->dam_params_size = 18;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 6;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 23:	num_global_params = 4;
        globals->uses_dam = 1;
        globals->num_params = 10;		//Need 1 extra for orifice_area
        globals->dam_params_size = 18;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 6;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 30:	num_global_params = 7;
        globals->uses_dam = 0;
        globals->num_params = 21;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 40:	num_global_params = 4;
        globals->uses_dam = 1;
        globals->num_params = 9;
        globals->dam_params_size = 6;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 6;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 60:	num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 4;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 101:	num_global_params = 0;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 105:	num_global_params = 0;
        globals->uses_dam = 0;
        globals->num_params = 16;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 190:	num_global_params = 6;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;
        //--------------------------------------------------------------------------------------------
    case 191:
        num_global_params = 7;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 3;
        break;
        //--------------------------------------------------------------------------------------------
    case 192:
        num_global_params = 7;
        globals->uses_dam = 0;
        globals->num_params = 7;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 3;
        break;
		//--------------------------------------------------------------------------------------------
	case 193:
		num_global_params = 0;
		globals->uses_dam = 0;
		globals->num_params = 4;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;
	break;
	//--------------------------------------------------------------------------------------------
	case 194:
		num_global_params = 3;
		globals->uses_dam = 0;
		globals->num_params = 4;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;
	break;
	//--------------------------------------------------------------------------------------------

	case 195:
		num_global_params = 5;
        globals->uses_dam = 0;
        globals->num_params = 6;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 3;
	break;
    //--------------------------------------------------------------------------------------------
    case 196:	num_global_params = 5;
        globals->uses_dam = 0;
        globals->num_params = 6;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 4;
        globals->min_error_tolerances = 5;
        break;
    //--------------------------------------------------------------------------------------------
    case 200:	num_global_params = 10;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 2;
        globals->areah_idx = 1;
        globals->num_disk_params = 12;
        globals->convertarea_flag = 1;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 219:	num_global_params = 3;
        globals->uses_dam = 0;
        globals->num_params = 5;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;
        //--------------------------------------------------------------------------------------------
    case 225:	num_global_params = 8;
        globals->uses_dam = 0;
        globals->num_params = 11;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 5;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;

    
    case 601:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 14;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 15;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;

    case 602:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 16;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 17;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break; 
    
    case 603:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 20;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 21;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;

    
    case 604:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 15;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 15;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;

    case 605:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 22;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 19;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;
 
    case 606:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 18;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 18;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;

    case 608:	num_global_params = 4;
        globals->uses_dam = 0;
        globals->num_params = 18;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 18;
        globals->convertarea_flag = 0;
        globals->num_forcings = 5;
        globals->min_error_tolerances = 5;
        break;
    case 609:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 18;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 18;
        globals->convertarea_flag = 0;
        globals->num_forcings = 4;
        globals->min_error_tolerances = 4;
        break;
    case 610:	num_global_params = 4;
        globals->uses_dam = 0;
        globals->num_params = 18;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 18;
        globals->convertarea_flag = 0;
        globals->num_forcings = 5;
        globals->min_error_tolerances = 4;
        break;
    case 607:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 15;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 15;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;

    case 654:	num_global_params = 1;
        globals->uses_dam = 0;
        globals->num_params = 14;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 14;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 3;
        break;
//--------------------------------------------------------------------------------------------
    case 249:	num_global_params = 12;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 7;
        break;
    case 250:	num_global_params = 9;
        globals->uses_dam = 0;
        globals->num_params = 9;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 4;
        globals->convertarea_flag = 0;
        globals->num_forcings = 1;
        globals->min_error_tolerances = 1;	//This should probably be higher...
        break;

//--------------------------------------------------------------------------------------------
    case 251:	num_global_params = 11;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 5;
        globals->min_error_tolerances = 4;
        break;
        //--------------------------------------------------------------------------------------------
    case 252:	num_global_params = 11;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 4;
        break;
        //--------------------------------------------------------------------------------------------
    case 253:	num_global_params = 11;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 4;
        globals->min_error_tolerances = 4;
        break;
        //--------------------------------------------------------------------------------------------
    case 254:	num_global_params = 12;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 7;
        break;
        //--------------------------------------------------------------------------------------------
    case 255:	num_global_params = 3;
        globals->uses_dam = 1;
        globals->num_params = 16;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 11;
        globals->convertarea_flag = 0;
        globals->num_forcings = 4;
        globals->min_error_tolerances = 5;
        break;
        //--------------------------------------------------------------------------------------------
    case 256:   num_global_params = 13;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 8;
        break;
        //--------------------------------------------------------------------------------------------
    case 257:   num_global_params = 31;
        globals->uses_dam = 0;
        globals->num_params = 9;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 4;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 8;
        break;
        //--------------------------------------------------------------------------------------------
    case 258:	num_global_params = 9;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;              //
        globals->areah_idx = 2;             //
        globals->num_disk_params = 3;       //
        globals->convertarea_flag = 0;      //
        globals->num_forcings = 4;          // runoff, infilt, evapot, reserv
        globals->min_error_tolerances = 7;  //
        break;
        //--------------------------------------------------------------------------------------------
    case 259:	num_global_params = 10;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;              //
        globals->areah_idx = 2;             //
        globals->num_disk_params = 3;       //
        globals->convertarea_flag = 0;      //
        globals->num_forcings = 4;          // runoff, infilt, evapot, reserv
        globals->min_error_tolerances = 7;  //
        break;
        //--------------------------------------------------------------------------------------------
    case 260:	num_global_params = 11;
        globals->uses_dam = 0;
        globals->num_params = 6;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 2;
        globals->min_error_tolerances = 4;
        break;
        //--------------------------------------------------------------------------------------------
    case 261:	num_global_params = 6;
        globals->uses_dam = 1;
        globals->num_params = 13;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 9;
        globals->convertarea_flag = 0;
        globals->num_forcings = 4;
        globals->min_error_tolerances = 8;
        break;
        //--------------------------------------------------------------------------------------------
    case 262:	num_global_params = 6;
        globals->uses_dam = 1;
        globals->num_params = 14;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 10;
        globals->convertarea_flag = 0;
        globals->num_forcings = 3;
        globals->min_error_tolerances = 8;
        break;
        //--------------------------------------------------------------------------------------------
  case 263: //256 with distributed params
		num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 14;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 15;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 8;
		break;
		        //--------------------------------------------------------------------------------------------
	case 264:   num_global_params = 13;
        globals->uses_dam = 0;
        globals->num_params = 8;
        globals->dam_params_size = 0;
        globals->area_idx = 0;
        globals->areah_idx = 2;
        globals->num_disk_params = 3;
        globals->convertarea_flag = 0;
        globals->num_forcings = 4;
        globals->min_error_tolerances = 8;
        break;
		//--------------------------------------------------------------------------------------------
	case 300:
		num_global_params = 6;
		globals->uses_dam = 0;
		globals->num_params = 20;
		globals->dam_params_size = 0;
		globals->area_idx = 2;
		globals->areah_idx = 1;
		globals->num_disk_params = 12;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;	//This should probably be higher...
		break;
		//--------------------------------------------------------------------------------------------
	case 301:
		num_global_params = 6;
		globals->uses_dam = 0;
		globals->num_params = 20;
		globals->dam_params_size = 0;
		globals->area_idx = 2;
		globals->areah_idx = 1;
		globals->num_disk_params = 12;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;	//This should probably be higher...
		break;
		//--------------------------------------------------------------------------------------------
	case 315:
		num_global_params = 6;
		globals->uses_dam = 0;
		globals->num_params = 20;
		globals->dam_params_size = 0;
		globals->area_idx = 2;
		globals->areah_idx = 1;
		globals->num_disk_params = 12;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;	//This should probably be higher...
		break;
		//--------------------------------------------------------------------------------------------
	case 400://tetis01
		num_global_params = 11;//v0,l1,l2,hu,infil,perc,surfvel,subrestime,gwrestime,meltfactor,tempthres
		globals->uses_dam = 0;
		globals->num_params = 6;//ai,li,ah,invtau,c1,c2
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 5; //precip, et, temperature,soil temperature,discharge
		globals->min_error_tolerances = 6; //as many as states:static,surface,subsurf,gw,channel,snow,
		break;
	case 612://tetis01
		num_global_params = 13;
		globals->uses_dam = 0;
		globals->num_params = 7;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 5;
		globals->min_error_tolerances = 7; //as many as states
		break;
    	//--------------------------------------------------------------------------------------------
	case 401://tetis02
		num_global_params = 11;
		globals->uses_dam = 0;
		globals->num_params = 6;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 5; //precip, et, temperature,soil temperature,discharge
		globals->min_error_tolerances = 10;  //as many as states
		break;
		//--------------------------------------------------------------------------------------------
	case 402://tetis03
		num_global_params = 4;
		globals->uses_dam = 1;
		globals->num_params = 11;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 3;
		globals->min_error_tolerances = 8;
		break;
		//--------------------------------------------------------------------------------------------
	case 2000:
		num_global_params = 6;
		globals->uses_dam = 0;
		globals->num_params = 20;
		globals->dam_params_size = 0;
		globals->area_idx = 2;
		globals->areah_idx = 1;
		globals->num_disk_params = 12;
		globals->convertarea_flag = 1;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;	//This should probably be higher...
		break;
		//--------------------------------------------------------------------------------------------
	default:
		printf("Error: Invalid model_uid (%u) in SetParamSizes.\n", model_uid);
		MPI_Abort(MPI_COMM_WORLD, 1);
		//--------------------------------------------------------------------------------------------
	}
    
  //Make sure the appropriate number of global parameters are given
	if (globals->num_global_params < num_global_params) {
		printf(
				"\nError: Obtained %u parameters from .gbl file. Expected %u for model model_uid %hu.\n",
				globals->num_global_params, num_global_params, model_uid);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (globals->num_global_params > num_global_params)
		printf(
				"\nWarning: Obtained %u parameters from .gbl file. Expected %u for model model_uid %hu.\n",
				globals->num_global_params, num_global_params, model_uid);
}

//Sets the function to be used when writing outputs. This method should set the following field:
//output_constrains_hdf5
//output_constrains_psql
//output_constrains_rec
void SetOutputConstraints(GlobalVars* globals)
{
    unsigned short int model_uid = globals->model_uid;
    //Set dim and start of differential variables
    switch (model_uid)
    {
        //--------------------------------------------------------------------------------------------
        case 196:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model196_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
        case 249:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model249_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;    
        case 254:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model254_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
        case 256:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model256_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
        case 264:
            globals->OutputConstrainsHdf5 = &OutputConstraints_Model256_Hdf5;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
        default:
            globals->OutputConstrainsHdf5 = NULL;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
    }
}

//Performs some unit conversions on the data in params. This takes place immediately after reading in the DEM data,
//so these changes are available in all routines of definetype.c, if params is available. Note that dam data
//and precalculations are not available here.
//VEC* params:		Vector of parameters to convert.
//unsigned int model_uid:	The index of the model.
void ConvertParams(
    double *params,
    unsigned int model_uid,
    void* external)
{
    if (model_uid == 19)
    {
        params[1] *= 1000;	//L: km -> m
        params[2] *= 1e6;	//A_h: km^2 -> m^2
    }
    else if (model_uid == 190 || model_uid == 191 || model_uid == 192 ||  model_uid==193 || model_uid==194 || model_uid == 195 || model_uid == 196)
    {
        params[1] *= 1000;	//L: km -> m
        params[2] *= 1e6;	//A_h: km^2 -> m^2
    }
    else if (model_uid == 20)
    {
        //params[0] *= 1e6;	//km^2 -> m^2
        params[1] *= 1000;	//km -> m
        params[2] *= 1e6;	//km^2 -> m^2
    }
    else if (model_uid == 60)
    {
        params[1] *= 1000;	//L: km -> m
        params[2] *= 1e6;	//A_h: km^2 -> m^2
        params[3] *= 1e-3;	//h_b: mm->m
    }
    else if (model_uid == 21)
    {
        //params[0] *= 1e6;	//km^2 -> m^2
        params[1] *= 1000;	//km -> m
        params[2] *= 1e6;	//km^2 -> m^2
    }
    else if (model_uid == 22 || model_uid == 23 || model_uid == 40)
    {
        //params[0] *= 1e6;	//km^2 -> m^2
        params[1] *= 1000;	//km -> m
        params[2] *= 1e6;	//km^2 -> m^2
    }
    else if (model_uid <= 5)
    {
        params[0] *= 1000;	//km -> m
        params[3] *= .001;	//mm -> m
        params[4] *= .001;	//mm -> m
    }
    else if (model_uid == 6)
    {
        params[0] *= 1000;	//km -> m
        params[3] *= .001;	//mm -> m
    }
    else if (model_uid == 15 || model_uid == 315)
    {
        params[0] *= 1000;	//L: km -> m
        params[3] *= .001;	//h_b: mm -> m
        params[4] *= .001;	//h_H: mm -> m
    }
    else if (model_uid == 30)
    {
        params[0] *= 1000;		//L_h:  km -> m
        params[4] *= .001;		//H_h:  mm -> m
        params[5] *= 1000.0; 	//MaxInfRate:  m/hr -> mm/hr
    }
    else if (model_uid == 105)
    {
        params[0] *= 1000;	//km -> m
        params[3] *= .001;	//mm -> m
        params[4] *= .001;	//mm -> m
    }
    else if (model_uid == 200)	//!!!! Did I screw these up on accident?!!!!
    {
        params[0] *= 1000;	//L_h:  km -> m
                                /*
                                //params[3] *= .001;	//mm -> m
                                params[4] *= .001;	//H_h:  mm -> m
                                params[5] *= 1000.0; //MaxInfRate:  m/hr -> mm/hr
                                */
    }
    else if (model_uid == 219 || model_uid == 225)
    {
        params[1] *= 1000;	//L: km -> m
        params[2] *= 1e6;	//A_h: km^2 -> m^2
    }
    else if (model_uid == 601 || model_uid == 602 || model_uid == 603 || model_uid == 604 || model_uid == 605 || model_uid == 606 || model_uid == 607)
    {
        params[1] *= 1000; //L: km -> m
        params[2] *= 1e6; // Ah: km^2 -> m^2
    }
    else if (model_uid == 654 || model_uid == 608 || model_uid == 609 || model_uid == 610)
    {
        params[1] *= 1000; //L: km -> m
        params[2] *= 1e6; // Ah: km^2 -> m^2
    } 
    else if (model_uid == 250)
    {
        params[1] *= 1000;		//L_h: km -> m
        params[2] *= 1e6;		//A_h: km^2 -> m^2
        params[4] *= .001;		//H_h: mm -> m
    }
    else if (model_uid == 249 || model_uid == 251 || model_uid == 252 || model_uid == 253 || 
    model_uid == 254 || model_uid == 255 || model_uid == 256 || model_uid == 257 || model_uid == 258 ||
     model_uid == 259 || model_uid == 260 || model_uid == 261 || model_uid == 262 || model_uid == 263 || 
     model_uid == 264 || model_uid==400 || model_uid==401 || model_uid == 612)
    {
        params[1] *= 1000;		//L_h: km -> m
        params[2] *= 1e6;		//A_h: km^2 -> m^2
    }
    else if (model_uid == 300 || model_uid == 301)
    {
        params[0] *= 1000;	//km -> m
        params[3] *= .001;	//mm -> m
        params[4] *= .001;	//mm -> m
    }
    else if (model_uid == 2000)
    {
        params[0] *= 1000;	//km -> m
        params[3] *= .001;	//mm -> m
        params[4] *= .001;	//mm -> m
    }
}

//Sets the system of ODEs and the Runge-Kutta solver for link. This method MUST set both link->differential
//	and link->solver. The Jacobian of f (link->jacobian) may be set here, if using an
//	implicit solver.
//Link* link: 		The link at which the ODEs and Runge-Kutta solver are selected.
//unsigned int model_uid: 	The index of the model to be set.
//unsigned int exp_imp: 0 if using an explicit solver, 1 if implicit.
//unsigned int dam: 	0 if no dam is present at link, 1 if a dam is present.
void InitRoutines(
    Link* link,
    unsigned int model_uid,
    unsigned int exp_imp,
    unsigned short dam,
    void* external)
{
    //Select appropriate RK Solver for the numerical method (link->solver)
    if ((model_uid == 21 || model_uid == 22 || model_uid == 23 || model_uid == 40 || model_uid == 261 || model_uid == 262) && dam == 1)
        link->solver = &ExplicitRKIndex1SolverDam;
    else if ((model_uid == 21 || model_uid == 22 || model_uid == 23 || model_uid == 40 || model_uid == 261 || model_uid == 262) && dam == 0)
        link->solver = &ExplicitRKIndex1Solver;
    else if (exp_imp == 0)
        link->solver = &ExplicitRKSolver;
    //	else if(link->method->exp_imp == 1)
    //		link->solver = &RadauRKSolver;
    else
        printf("Warning: No solver selected for link ID %u.\n", link->ID);

    //Select the RHS function of the ODE (link->differential, link->jacobian)
    if (model_uid == 0)
    {
        link->dim = 1;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &simple_river;
        link->jacobian = &Jsimple;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_1States;
    }
    else if (model_uid == 1 || model_uid == 3)
    {
        link->dim = 2;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &river_rainfall;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_2States;
    }
    else if (model_uid == 2)
    {
        link->dim = 2;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &simple_hillslope;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_2States;
    }
    else if (model_uid == 4)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &simple_soil;
        link->jacobian = &Jsimple_soil;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 5)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &soil_rainfall;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Model5;
    }
    else if (model_uid == 6)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &qsav_rainfall;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 15)
    {
        link->dim = 2;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &river_rainfall_adjusted;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_2States;
    }
    else if (model_uid == 19)
    {
        link->dim = 3;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &LinearHillslope_Evap;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_3States;
    }
    else if (model_uid == 20)
    {
        link->dim = 3;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &Hillslope_Toy;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_3States;
    }
    else if (model_uid == 21)
    {
        link->dim = 4;
        link->no_ini_start = 2;
        link->diff_start = 1;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 1;

        if (dam)	link->differential = &dam_rain_hillslope;
        else	link->differential = &nodam_rain_hillslope;
        link->algebraic = &dam_q;
        link->check_state = &dam_check;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 22)
    {
        link->dim = 4;
        link->no_ini_start = 2;
        link->diff_start = 1;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 1;

        if (dam)	link->differential = &dam_rain_hillslope2;
        else	link->differential = &nodam_rain_hillslope2;
        link->algebraic = &dam_q2;
        link->check_state = &dam_check2;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 23)
    {
        link->dim = 4;
        link->no_ini_start = 2;
        link->diff_start = 1;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 1;

        if (dam)	link->differential = &dam_rain_hillslope3;
        else	link->differential = &nodam_rain_hillslope3;
        link->algebraic = &dam_q3;
        link->check_state = &dam_check3;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 30)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &lcuencas_soilrain;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Model30;
    }
    else if (model_uid == 40)
    {
        link->dim = 4;
        link->no_ini_start = 2;
        link->diff_start = 1;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 1;

        if (dam)	link->differential = &dam_rain_hillslope_qsv;
        else	link->differential = &nodam_rain_hillslope_qsv;
        link->algebraic = &dam_q_qvs;
        link->check_state = &dam_check_qvs;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 60)
    {
        link->dim = 3;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &LinearHillslope_Evap_RC;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_3States;
    }
    else if (model_uid == 101)
    {
        link->dim = 3;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &Robertson;
        link->jacobian = &JRobertson;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_3States;
    }
    else if (model_uid == 105)
    {
        link->dim = 2;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &river_rainfall_summary;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_2States;
    }
    else if (model_uid == 190)
    {
        link->dim = 3;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &LinearHillslope_MonthlyEvap;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_3States;
    }
    else if (model_uid == 191)
    {
        link->dim = 6;
        link->no_ini_start = 3;
        link->diff_start = 0;

        link->num_dense = 2;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 5;

        if (link->has_res)
        {
            link->differential = &LinearHillslope_Reservoirs_extras;
            link->solver = &ForcedSolutionSolver;
        }
        else	link->differential = &LinearHillslope_MonthlyEvap_extras;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 192)
    {
        link->dim = 6;
        link->no_ini_start = 3;
        link->diff_start = 0;

        link->num_dense = 2;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 5;

        if (link->has_res)
        {
            link->differential = &LinearHillslope_Reservoirs_extras;
            link->solver = &ForcedSolutionSolver;
        }
        else	link->differential = &LinearHillslope_MonthlyEvap_kI_extras;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
        else if (model_uid == 193)
    {
        link->dim = 1;
        link->no_ini_start = 1;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &routing_runoff2;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 194)
    {
        link->dim = 1;
        link->no_ini_start = 1;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &routing_runoff1;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
	else if (model_uid == 195)
    {
        link->dim = 4;
        link->no_ini_start = 3;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &LinearHillslope_MonthlyEvap_OnlyRouts;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
	else if (model_uid == 196)
	{
		link->dim = 5;
		link->no_ini_start = 3;
		link->diff_start = 0;

		link->num_dense = 1;

		link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;

		if (link->has_res)
		{
			link->differential = &LinearHillslope_MonthlyEvap_OnlyRouts_HasReservoir;
			link->solver = &ForcedSolutionSolver;
		}
		else
		{
			link->differential = &LinearHillslope_MonthlyEvap_OnlyRouts_NotReservoir;
		}
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
    else if (model_uid == 200)	//This is for use with SIMPLE only
    {
        link->dim = 2;
        link->no_ini_start = 1;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = NULL;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_2States;
    }
    else if (model_uid == 219)
    {
        link->dim = 1;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &Tiling;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_1States;
    }
    else if (model_uid == 225)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &TilingHillslope;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    
    else if (model_uid == 601)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->differential = &VariableHillslope;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    

    else if (model_uid == 602)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->differential = &VariableInterflow;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }


    else if (model_uid == 603)
    {
        link->dim = 5;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->differential = &VariableTriLayer;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_5States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }


    else if (model_uid == 604)
    {
        link->dim = 7;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->differential = &VariableThreshold;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
            //&CheckConsistency_Nonzero_5States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }


    else if (model_uid == 605)
    {
        link->dim = 7;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->differential = &VariableThreshold3;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
            //&CheckConsistency_Nonzero_5States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }

    
    else if (model_uid == 606)
    {
        link->dim = 7;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->differential = &ExponentialExp;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
            //&CheckConsistency_Nonzero_5States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }

    else if (model_uid == 608)
    {
        link->dim = 5;
        link->no_ini_start = 5; //link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        //link->dense_indices[1] = 4;
        //link->dense_indices[2] = 5;

        if (link->has_res)
        { 
            link->differential = &Tiles_Reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else    link->differential = &TilesModel;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    
    else if (model_uid == 609)
    {
        link->dim = 5;
        link->no_ini_start = 5; //link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        if (link->has_res)
        { 
            link->differential = &Tiles_Reservoirs_Base;
            link->solver = &ForcedSolutionSolver;
        }
        else    link->differential = &TilesModel_Base;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }

    else if (model_uid == 610)
    {
        link->dim = 5;
        link->no_ini_start = 5; //link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        if (link->has_res)
        { 
            link->differential = &ActiveLayerSnow_Reservoir;
            link->solver = &ForcedSolutionSolver;
        }
        else    link->differential = &ActiveLayerSnow;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }

    else if (model_uid == 654)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->differential = &Variable_TopLayer;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
        //link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 249)
    {
        link->dim = 6;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 3;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 4;
        link->dense_indices[2] = 5;

        if (link->has_res)
        {
            link->differential = &model249_reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else			
            link->differential = &model249;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 250)
    {
        link->dim = 3;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &NonLinearHillslope;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_3States;
    }
            else if (model_uid == 251)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        if (link->has_res)
        {
            link->differential = &TopLayerHillslope_Reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else			link->differential = &model251;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 252)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &model252;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 253)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        if (link->has_res)
        {
            link->differential = &TopLayerHillslope_Reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else			link->differential = &model253;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 254)
    {
        link->dim = 7;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 2;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 6;

        if (link->has_res)
        {
            link->differential = &TopLayerHillslope_Reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else			
            link->differential = &model254;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 255)
    {
        link->dim = 5;
        link->no_ini_start = 5;
        link->diff_start = 1;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));

        if (link->has_res)
        {
            link->dense_indices[0] = 0;	//Discharge
            link->differential = &TopLayerHillslope_Reservoirs; //common to various models
            link->solver = &ForcedSolutionSolver;
        }
        else
        {
            link->dense_indices[0] = 1;	//Storage
            link->differential = &model255;
            link->solver = &ExplicitRKIndex1SolverDam;
        }
        link->algebraic = &dam_model255;
        link->check_state = &dam_check_qvs;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 256)
    {
        link->dim = 9;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 3;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 7;
        link->dense_indices[2] = 8;

        if (link->has_res)
        {
            link->differential = &model256_reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else		
            link->differential = &model256;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
	    
    else if (model_uid == 257)
    {
        link->dim = 8;
        link->no_ini_start = 4;
        link->diff_start = 0;

        link->num_dense = 2;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 7;

        if (link->has_res)
        {
            link->differential = &TopLayerHillslope_Reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else
            link->differential = &TopLayerHillslope_spatial_velocity;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if ((model_uid == 258)||(model_uid == 259))
    {
        link->dim = 8;
        link->no_ini_start = 4;
        link->diff_start = 0;

        link->num_dense = 2;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 6;

        if(model_uid == 258)
        {
            if (link->has_res)
            {
                link->differential = &OfflineTopLayerHillslope_Reservoirs;
                link->solver = &ForcedSolutionSolver;
            }
            else
            {
                link->differential = &OfflineTopLayerHillslope;
            }
        }
        else if (model_uid == 259)
        {
            if (link->has_res)
            {
                link->differential = &OfflineTopLayerInterflowHillslope_Reservoirs;
                link->solver = &ForcedSolutionSolver;
            }
            else
            {
                link->differential = &OfflineTopLayerInterflowHillslope;
            }
        }
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 260)
    {
        link->dim = 4;
        link->no_ini_start = link->dim;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &TopLayerNonlinearExp;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_4States;
    }
    else if (model_uid == 261)
    {
        link->dim = 8;
        link->no_ini_start = 5;
        link->diff_start = 1;

        link->num_dense = 2;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 1;
        link->dense_indices[1] = 7;

        if (link->has_res)
        {
            link->differential = &TopLayerNonlinearExpSoilvel_Reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else
        {
            link->differential = &TopLayerNonlinearExpSoilvel;
        }
        link->algebraic = &dam_TopLayerNonlinearExpSoilvel;
        link->check_state = &dam_check_qvs;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    else if (model_uid == 262)
    {
        link->dim = 8;
        link->no_ini_start = 5;
        link->diff_start = 1;

        if (link->has_res)
        {
            link->num_dense = 1;
            link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
            link->dense_indices[0] = 0;

            link->differential = &TopLayerNonlinearExpSoilvel_ConstEta_Reservoirs;
            link->algebraic = NULL;
            link->solver = &ForcedSolutionSolver;
        }
        else
        {
            link->num_dense = 2;
            link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
            link->dense_indices[0] = 1;
            link->dense_indices[1] = 7;
            link->differential = &TopLayerNonlinearExpSoilvel_ConstEta;
            link->algebraic = &dam_TopLayerNonlinearExpSoilvel_ConstEta;
        }
        link->check_state = &dam_check_qvs;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_qs;
    
	} else if (model_uid == 263) {
		link->dim = 8;
		link->no_ini_start = 4;
		link->diff_start = 0;

		link->num_dense = 2;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices,
				link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		link->dense_indices[1] = 7;

		if (link->has_res) {
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		} else
			link->differential = &model263;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
    else if (model_uid == 264)
    {
        link->dim = 8;
        link->no_ini_start = 4;
        link->diff_start = 0;

        link->num_dense = 2;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;
        link->dense_indices[1] = 7;

        if (link->has_res)
        {
            link->differential = &TopLayerHillslope_Reservoirs;
            link->solver = &ForcedSolutionSolver;
        }
        else		
            link->differential = &model264;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
	else if (model_uid == 612) //tetis01
			{
		link->dim = 7;
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices,
				link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		//link->dense_indices[1] = 7;

		if (link->has_res) {
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		} else
			link->differential = &tetis_nicoV1;
		    link->algebraic = NULL;
		    link->check_state = NULL;
		    link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	} 
    else if (model_uid == 400) //tetis01
			{
		link->dim = 6;
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices,
				link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		//link->dense_indices[1] = 7;

		if (link->has_res) {
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		} else
			link->differential = &model400;
		    link->algebraic = NULL;
		    link->check_state = NULL;
		    link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	} 
    else if (model_uid == 401) //tetis02
	{
		link->dim = 10;//q,basinrain,basinsurf,basinsubsur,basingw,static,surf,subsurf,gw,snow
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices,
		link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		//link->dense_indices[1] = 7;

		if (link->has_res) {
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		} else
			link->differential = &model401;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	} 
	else if (model_uid == 402) //tetis03
				{
			link->dim = 9;
			link->no_ini_start = 4;
			link->diff_start = 0;

			link->num_dense = 2;
			link->dense_indices = (unsigned int*) realloc(link->dense_indices,
					link->num_dense * sizeof(unsigned int));
			link->dense_indices[0] = 0;
			link->dense_indices[1] = 7;

			if (link->has_res) {
				link->differential = &Tetis03_Reservoirs;
				link->solver = &ForcedSolutionSolver;
			} else
				link->differential = &Tetis03;
			link->algebraic = dam_Tetis03;
			link->check_state = NULL;
			link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
		}
	//else if (model_uid == 300)
	//{
	//    link->dim = 2;
	//    link->no_ini_start = 1;
	//    link->diff_start = 0;

	//    link->num_dense = 1;
	//    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
	//    link->dense_indices[0] = 0;

	//    link->differential = &assim_simple_river;
	//    link->algebraic = NULL;
	//    link->check_state = NULL;
	//    link->check_consistency = &CheckConsistency_Nonzero_1States;
	//}
	//else if (model_uid == 301)	//!!!! For data assimilation. Needs updating. !!!!
	//{
	//    printf("!!!! InitRoutines: model 301 needs to be updated. !!!!\n");

	//    link->dim = 3;
	//    link->no_ini_start = 2;
	//    link->diff_start = 0;

	//    link->num_dense = 1;
	//    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
	//    link->dense_indices[0] = 0;

	//    link->differential = &assim_river_rainfall;
	//    link->algebraic = NULL;
	//    link->check_state = NULL;
	//    link->check_consistency = &CheckConsistency_Nonzero_2States;
	//}
	//else if (model_uid == 315)
	//{
	//    printf("!!!! InitRoutines: model 315 needs to be updated. !!!!\n");

	//    link->dim = 2;
	//    link->no_ini_start = 2;
	//    link->diff_start = 0;

	//    link->num_dense = 1;
	//    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
	//    link->dense_indices[0] = 0;

	//    link->differential = &assim_river_rainfall_adjusted;
	//    link->algebraic = NULL;
	//    link->check_state = NULL;
	//    link->check_consistency = &CheckConsistency_Nonzero_2States;
	//}
	/*
	 else if(model_uid == 2000)
	 {
	 link->differential = &parser_test;
	 link->algebraic = NULL;
	 link->check_state = NULL;
	 }
	 */
	else
		printf("Warning: No ODE selected for link ID %u.\n", link->ID);
}

//Perform precalculations needed for the differential equation.  These should be stored in params after the DEM
//	data and after the dam data (i.e. params[disk_params] is the first precalcuation, params[params_size]
//	is the first dam datum). This method is run at each link. This method can do nothing, if no precalculations are
//	expected for the model. This assumes the same number of precalculations regardless if there is a dam or not.
//VEC* global_params:		Vector of the global parameters of the system. These are already set and are available for use.
//VEC* params:			Vector of the parameters at a link. The DEM data and dam data are already set (but may be
//					altered here). Only the entries for precalculations need to be set.
//unsigned int disk_params:	The first entry of params that should be set here.
//unsigned int params_size:	First entry of the dam data. Don't change this entry or later unless you want to modify the dam!
//unsigned int model_uid:		The index of the model.

void Precalculations(
    Link* link_i,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_disk_params, unsigned int num_params,
    unsigned short dam,
    unsigned int model_uid,
    void* external)
{
    if (model_uid == 19)
    {
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
        //The numbering is:	0   1   2   3  4    5    6   7
        //Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g,e_pot
        //The numbering is:        0      1        2     3  4   5   6
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double RC = global_params[3];
        double v_h = global_params[4];
        double v_g = global_params[5];

        vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
        vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
        vals[5] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[6] = RC*(0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[7] = (1.0 - RC)*(0.001 / 60.0);	//(mm/hr->m/min)  c_2
    }
    else if ((model_uid == 190) || (model_uid == 191))
    {
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
        //The numbering is:	    0   1   2   3  4    5    6   7
        //Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g (,v_B)
        //The numbering is:        0      1        2     3  4   5     6
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double RC = global_params[3];
        double v_h = global_params[4];
        double v_g = global_params[5];

        vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
        vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
        vals[5] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[6] = RC*(0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[7] = (1.0 - RC)*(0.001 / 60.0);	//(mm/hr->m/min)  c_2
    }
    else if (model_uid == 192)
    {
        //Order of parameters: A_i,L_i,A_h,k_2,k_i,invtau,c_1
        //The numbering is:	    0   1   2   3   4    5     6
        //Order of global_params: v_r,lambda_1,lambda_2,k_I_factor,v_h,k_3,v_B
        //The numbering is:        0      1        2        3       4   5   6
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double k_i_factor = global_params[3];
        double v_h = global_params[4];

        vals[3] = v_h * L_i / A_h * 60.0;                                // [1/min]  k_2
        vals[4] = vals[3] * k_i_factor;	                                 // [1/min]  k_i
        vals[5] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);  // [1/min]  invtau
        vals[6] = 0.001 / 60.0;                                          // (mm/hr->m/min)  c_1
    }
    else if (model_uid == 193)
    {
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
        //Order of parameters: A_i,L_i,A_h,invtau
        //The numbering is:	0   1   2   3
        //Order of global_params: v_r,lambda_1,lambda_2,
        //The numbering is:        0      1        2
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        //double v_r = global_params[0];
        //double lambda_1 = global_params[1];
        //double lambda_2 = global_params[2];
        //double v_h = global_params[3];
        //double k3 = global_params[4];

        //vals[3] = v_h * L_i / A_h * 60.0;	                            // [1/min]  k2
        //vals[4] = k3;                                                   // [1/min]  k3
        //vals[3] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }
    else if (model_uid == 194)
    {
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
        //Order of parameters: A_i,L_i,A_h,invtau
        //The numbering is:	0   1   2   3
        //Order of global_params: v_r,lambda_1,lambda_2,
        //The numbering is:        0      1        2
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        //double v_h = global_params[3];
        //double k3 = global_params[4];

        //vals[3] = v_h * L_i / A_h * 60.0;	                            // [1/min]  k2
        //vals[4] = k3;                                                   // [1/min]  k3
        vals[3] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }
	else if (model_uid == 195)
    {
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
        //The numbering is:	0   1   2   3  4    5    6   7
        //Order of global_params: v_r,lambda_1,lambda_2,v_h,k_3 (,v_B)
        //The numbering is:        0      1        2     3   4     5
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double k3 = global_params[4];

        vals[3] = v_h * L_i / A_h * 60.0;	                            // [1/min]  k2
        vals[4] = k3;                                                   // [1/min]  k3
        vals[5] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }
    else if (model_uid == 196)
    {
        //Order of parameters: A_i,L_i,A_h,k_2,k_3,invtau,c_1,c_2
        //The numbering is:	    0   1   2   3   4     5    6   7
        //Order of global_params: v_r,lambda_1,lambda_2,v_h, k3 (,v_B)
        //The numbering is:        0      1        2     3   4     5
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double k3 = global_params[4];

        vals[3] = v_h * L_i / A_h * 60.0;	                            // [1/min]  k2
        vals[4] = k3;                                                   // [1/min]  k3
        vals[5] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }
    else if (model_uid == 20)
    {
        //Order of parameters: A_i,L_i,A_h,invtau,c_1,c_2
        //The numbering is:	0   1   2    3     4   5
        //Order of global_params: v_r,lambda_1,lambda_2,beta,k_p,k_a,theta_p,theta_a,scale_p,scale_a
        //The numbering is:        0      1        2     3    4   5     6       7      8        9
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];


        vals[3] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i); //[1/min]  invtau

        /*
        //For gamma
        double k_p = global_params[4];
        double k_a = global_params[5];
        double theta_p = global_params[6];
        double theta_a = global_params[7];
        vals[4] = 1.0/(tgamma(k_p)*pow(theta_p,k_p)); //c_1
        vals[5] = 1.0/(tgamma(k_a)*pow(theta_a,k_a)); //c_2
        */
        /*
        //For log normal
        double mu = global_params[3];
        double sigma2 = global_params[4];
        double scale = global_params[5];
        vals[4] = scale/(pow(2.0*3.141592653589*sigma2,0.5));
        vals[5] = 0.0;
        */

        //Set some random values
        //vals[4] = (double)(rand()%996) / 10.0 + 0.5;	//0.5 to 100.0
        //vals[5] = (double)(rand()%991) + 10.0;	//10 to 1000

        /*
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
        //The numbering is:	0   1   2   3  4    5    6   7
        //Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
        //The numbering is:        0      1        2     3  4   5
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double RC = global_params[3];
        double v_h = global_params[4];
        double v_g = global_params[5];
        double F_et = 0.05;

        vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
        vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
        vals[5] = 60.0*v_r*pow(A_i,lambda_2) / ((1.0-lambda_1)*L_i);	//[1/min]  invtau
        vals[6] = RC*(0.001/60.0);		//(mm/hr->m/min)  c_1
        vals[7] = F_et*(1.0-RC)*(0.001/60.0);	//(mm/hr->m/min)  c_2
        */
    }
    else if (model_uid == 60)
    {
        //Order of parameters: A_i,L_i,A_h,h_b,k2,k3,invtau,c_1
        //The numbering is:	0   1   2   3   4  5    6    7 
        //Order of global_params: v_r,lambda_1,lambda_2,v_h,v_g,e_pot
        //The numbering is:        0      1        2     3   4   5  
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double v_g = global_params[4];

        vals[4] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
        vals[5] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
        vals[6] = 60.0*v_r*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
    }
    else if (model_uid == 21)
    {
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
        //The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
        //Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h,v_g
        //The numbering is:        0      1        2     3  4   5   6
        //Need to set entries 3, 4, 5, and 6 of params.
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[5];
        double v_g = global_params[6];

        vals[3] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
                                            //vals[4] = vals[3] / 20.0;		//[1/min]  k3
        vals[4] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
        vals[5] = pow(v_r*pow(A_i, lambda_2) / L_i, 1.0 / (1.0 - lambda_1)) * 60.0; //[1/min]  invtau

        if (dam)		vals[6] = 3.1415926535897932 * vals[11] * vals[11] / 4.0;
    }
    else if (model_uid == 22 || model_uid == 23)
    {
        //Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
        //The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
        //Order of global_params: lambda_1,lambda_2,S_0,v_g
        //The numbering is:         0        1       2   3
        //Need to set entries 6, 7, 8, and 9 of params.
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_h = params[4];
        double v_r = params[5];
        double lambda_1 = global_params[0];
        double lambda_2 = global_params[1];
        double v_g = global_params[3];

        vals[6] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
                                            //vals[7] = vals[6] / 20.0;		//[1/min]  k3
        vals[7] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
        vals[8] = pow(v_r*pow(A_i, lambda_2) / L_i, 1.0 / (1.0 - lambda_1)) * 60.0; //[1/min]  invtau

        if (dam)		vals[9] = 3.1415926535897932 * vals[14] * vals[14] / 4.0; //orifice_area
    }
    else if (model_uid == 40)
    {
        //Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
        //The numbering is:	0   1   2  3   4   5   6  7   8
        //Order of global_params: lambda_1,lambda_2,S_0,v_g
        //The numbering is:         0        1       2   3
        //Need to set entries 6, 7, and 8 of params.
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_h = params[4];
        double v_r = params[5];
        double lambda_1 = global_params[0];
        double lambda_2 = global_params[1];
        double v_g = global_params[3];

        vals[6] = v_h * L_i / A_h * 60.0;	//[1/min]  k2
        vals[7] = v_g * L_i / A_h * 60.0;	//[1/min]  k3
        vals[8] = pow(v_r*pow(A_i, lambda_2) / L_i, 1.0 / (1.0 - lambda_1)) * 60.0; //[1/min]  invtau
    }
    else if (model_uid <= 5 || model_uid == 200)
    {
        //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
        //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
        //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC, (for 200) u_0,U,DL,scale
        //The numbering is:        0      1        2     3   4   5   	       6  7  8   9
        //Need to set entries 12-19 of params.
        double* vals = params;
        double K_T = 1.0;
        double s_r = 1.0;
        double rootS_h = pow(vals[7], .5);
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double Q_r = global_params[3];
        double A_r = global_params[4];
        double RC = global_params[5];

        vals[12] = 60.0*v_r*pow(vals[2] / A_r, lambda_2) / ((1.0 - lambda_1)*vals[0]);
        vals[13] = vals[3] / s_r;
        vals[14] = 2.0 / .6 * vals[0] * s_r / Q_r * 1.0 / vals[8] * rootS_h * pow(vals[3], 2.0 / 3.0);	//c_1
        vals[15] = vals[6] * vals[0] * vals[3] / (3600.0 * Q_r);
        vals[16] = 1e-3 / (60 * s_r) * RC;	//c_3
        vals[17] = 60 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * 1.0 / vals[8] * rootS_h;	//c_4
        vals[18] = K_T / 60.0;
        vals[19] = vals[6] / (60.0*s_r);
    }
    else if (model_uid == 6)
    {
        //Order of parameters: L_i,A_h,A_i,h_b,K_sat,K_sp,d_0,d_1,d_2,invbeta,alpha_N,alpha_soil,a_res,v_res,invtau,gamma,c_1,c_2,c_3,c_4
        //The numbering is:     0   1   2   3   4     5    6   7   8     9      10       11       12    13     14    15   16  17  18  19
        //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r
        //The numbering is:        0      1        2     3   4
        //Need to set entries 14 - 19 and 9 of params.
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double Q_r = global_params[3];
        double A_r = global_params[4];

        double* vals = params;
        double L_i = vals[0];
        double A_h = vals[1];
        double A_i = vals[2];
        double h_b = vals[3];
        double K_sat = vals[4];
        double K_sp = vals[5];
        double beta = vals[9];
        double alpha_N = vals[10];
        double alpha_soil = vals[11];

        vals[9] = 1.0 / beta;		//invbeta

        vals[14] = 60.0 * v_r * pow(A_i / A_r, lambda_2) / ((1.0 - lambda_1) * L_i);	//invtau
        vals[15] = (1e6 / 60.0) * A_h / Q_r;	//gamma
        vals[16] = K_sp / (60.0 * h_b);	//c_1
        vals[17] = (1e-6 / 60.0) * alpha_soil * K_sat * L_i / (alpha_N * A_h); //c_2
        vals[18] = 1e-3 / 60.0; //c_3
        vals[19] = alpha_N * h_b; //c_4
    }
    else if (model_uid == 15)
    {
        //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
        //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
        //Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
        //The numbering is:        0      1        2     3   4   5
        //Need to set entries 12-19 of params.
        double* vals = params;
        double K_T = 1.0;
        double s_r = 1.0;
        double rootS_h = pow(vals[7], .5);
        double L = params[0];
        double A_h = params[1] * 1e6;	//Put into m^2
        double eta = params[8];
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double A_r = global_params[4];
        double RC = global_params[5];


        vals[12] = 60.0*v_r*pow(vals[2] / A_r, lambda_2) / ((1.0 - lambda_1)*vals[0]);	//invtau
        vals[13] = vals[3] / s_r; //epsilon
        vals[14] = v_h*L;	//c_1 [m^2/s]
        vals[15] = vals[6] * vals[0] * vals[3] / 3600.0; //c_2
        vals[16] = (1e-3 / 60.0) * RC;	//c_3
        vals[17] = 60.0*v_h*L / A_h;	//c_4 [1/min], A_h converted above
        vals[18] = K_T / 60.0;
        vals[19] = vals[6] / (60.0*s_r);
    }
    else if (model_uid == 30)
    {
        //y_i has q (y_i[0]), s_p (y_i[1]), h_w (y_i[2]), theta (y_i[3])
        //Order of parameters: L_h,A_h,A_up,H_b,H_h,max_inf_rate,K_SAT,S_H,n_vh, b, c, d,K_Q,V_T,c_1,c_2,c_3,c_4,c_5,c_6,c_7
        //The numbering is:     0   1   2    3   4       5         6    7   8    9 10 11  12  13  14  15  16  17  18  19  20
        //Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,K_T,C_r
        //The numbering is:        0      1        2     3   4   5   6
        //Need to set entries 12-19 of params.

        //Global params
        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double A_r = global_params[4];
        double C_r = global_params[6];

        //Local params
        double L_H = params[0];
        double A_H = params[1];
        double A_up = params[2];
        double H_b = params[3];
        double H_h = params[4];
        double K_SAT = params[6];
        double S_H = params[7];
        double n_vh = params[8];
        double H_relmax = H_h;


        params[12] = 60.0*C_r*v_0*pow(A_up / A_r, lambda_2) / ((1.0 - lambda_1)*L_H);	//K_Q
        params[13] = 1e3 * H_b * A_H;	//V_T
        params[14] = 3.6e3 / n_vh * pow(S_H, .5);	//c_1
                                                        //params[14] = 3.6e6 / n_vh * pow(S_H,.5);	//c_1
                                                        //params[15] = (2e-6/.6) * (L_H/A_H);	//c_2
        params[15] = (2e-6) * (L_H / A_H);	//c_2
        params[16] = 1e3 / params[13] * K_SAT * S_H;	//c_3
        params[17] = (1.0 / 3.6) * A_H;	//c_4
        params[18] = (1e6 / 60.0) * A_H / H_b;	//c_5
        params[19] = (1e3 / 60.0) * A_H;	//c_6
        params[20] = (H_relmax > 1e-12) ? 1e6 * A_H / H_relmax : 0.0;	//c_7
    }
    else if (model_uid == 105)
    {
        //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,c_1,c_2,c_3
        //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13  14  15 
        //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
        //The numbering is:        0      1        2     3   4   5
        double* vals = params;
        double v_r = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double A_r = global_params[4];
        double RC = global_params[5];
        double rootS_h = pow(vals[7], .5);

        //invtau
        vals[12] = 60.0*v_r*pow(vals[2] / A_r, lambda_2) / ((1.0 - lambda_1)*vals[0]);
        //c_1
        vals[13] = 2.0 / .6 * vals[0] * rootS_h / vals[8];
        //c_2
        vals[14] = 1e-3 / 60.0 * RC;
        //c_3
        vals[15] = 60.0 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * rootS_h / vals[8];
    }
    else if (model_uid == 219)
    {
        //Order of parameters: A_i,L_i,A_h,invtau,c_1
        //The numbering is:	0   1   2    3     4
        //Order of global_params: v_0,lambda_1,lambda_2
        //The numbering is:        0      1        2
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];

        vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[4] = (1e-3 / 3600.0) * A_h;		//c_1
    }
    else if (model_uid == 225)
    {
	//Order of parameters: A_i,L_i,A_h,slope,tileQ,invtau,aNo,actT,LinSl,aTile,LinIn,
	//The numbering is:	0   1   2    3     4     5     6   7     8   9     10
	//Order of global_params: v_r,lambda_1,lambda_2,tL,bL,kdry,ki,k2
	//The numbering is:        0      1        2     3 4   5   6   7 
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
	    double slope = params[3];
    	double tileQ = params[4];
        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];

        vals[5] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[6] = 1.1e-15 * slope + 1.4e-14;				// aNo
    	vals[7] = 0.005 * slope + 1.7;					// actTile
	    vals[8] = 2.7e-8 * slope + 0.0021;				// Linslope
    	vals[9] = (vals[6] * vals[7] * exp(10 * vals[7]))/(exp(21 * vals[7]));	// aTile
	    vals[10] = vals[6] * vals[7] * exp(10 * vals[7]) - vals[8] * vals[7];	// LinIntercept
    }
    else if (model_uid == 601)
    {
	//Order of parameters: A_i,L_i,A_h,slope,tileQ,invtau,aNo,actT,LinSl,aTile,LinIn,
	//The numbering is:	0   1   2    3     4     5     6   7     8   9     10
	//Order of global_params: v_r,lambda_1,lambda_2,tL,bL,kdry,ki,k2
	//The numbering is:        0      1        2     3 4   5   6   7 
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
	    double slope = params[3];
        double v_r = params[4];
        double a_r = params[5];
        double v_s = params[6];
        double a_s = params[7];
        double k1 = params[8];
        double k2 = params[9];
        double t_l = params[10];
        double b_L = params[11];
        double lambda_1 = params[12];
        double lambda_2 = params[13];
        double v_0 = params[14];   

        vals[13] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
    }
    else if (model_uid == 602)
    {
	//Order of parameters: A_i,L_i,A_h,slope,tileQ,invtau,aNo,actT,LinSl,aTile,LinIn,
	//The numbering is:	0   1   2    3     4     5     6   7     8   9     10
	//Order of global_params: v_r,lambda_1,lambda_2,tL,bL,kdry,ki,k2
	//The numbering is:        0      1        2     3 4   5   6   7 
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
	    double slope = params[3];
        double v_r = params[4];
        double a_r = params[5];
        double v_l = params[6];
        double a_l = params[7];
        double v_s = params[8];
        double a_s = params[9];
        double k1 = params[10];
        double k2 = params[11];
        double t_l = params[12];
        double b_L = params[13];
        double lambda_1 = params[14];
        double lambda_2 = params[15];
        double v_0 = params[16];   

        vals[15] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
    }
    else if (model_uid == 603)
    {
	//Order of parameters: A_i,L_i,A_h,slope,tileQ,invtau,aNo,actT,LinSl,aTile,LinIn,
	//The numbering is:	0   1   2    3     4     5     6   7     8   9     10
	//Order of global_params: v_r,lambda_1,lambda_2,tL,bL,kdry,ki,k2
	//The numbering is:        0      1        2     3 4   5   6   7 
        double* vals = params;
        double Ai = params[0];
        double Li = params[1];
        double Ah = params[2];
	    double slope = params[3];
        double v1 = params[4];
        double a1 = params[5];
        double v2 = params[6];
        double a2 = params[7];
        double v3 = params[8];
        double a3 = params[9]; 
        double v4 = params[10];
        double a4 = params[11];
        double h1 = params[12];
        double h2 = params[13];
        double h3 = params[14];
        double k1 = params[15];
        double k2 = params[16];
        double k3 = params[17];
        double lambda1 = params[18];
        double lambda2 = params[19];
        double v0 = params[20];   

        vals[19] = 60.0*v0*pow(Ai, lambda2) / ((1.0 - lambda1)*Li);	//[1/min]  invtau
    }

    else if (model_uid == 604)
    {
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = params[3];
        double a_r = params[4];
        double v_s1 = params[5];
        double v_s2 = params[6];
        double k1 = params[7];
        double k2 = params[8];
        double t_L = params[9];
        double Si = params[10];
        double Sn = params[11];
        double lambda_1 = params[12];
        double lambda_2 = params[13];
        double v0 = params[14];

        vals[13] = 60.0*v0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[14] = v_s1*(Si - Sn);//(v_s1 - v_s2)*b_L;  //second line intercept coefficient

    }

    else if (model_uid == 605)
    {
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = params[3];
        double a_r = params[4];
        double v_s1 = params[5];
        double v_s2 = params[6];
        double v_s3 = params[7];
        double v_s4 = params[8];
        double k1 = params[9];
        double k2 = params[10];
        double t_L = params[11];
        double S1 = params[12];
        double S2 = params[13];
        double S3 = params[14];
        double S4 = params[15];
        double lambda_1 = params[16];
        double lambda_2 = params[17];
        double v0 = params[18];

        vals[17] = 60.0*v0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[18] = v_r * (L_i / A_h) * 60; // [1/min] runoff speed.


        vals[19] = v_s1*(S2 - S1);// first intercept
        vals[20] = v_s2*(S3 - S2) + vals[19];// second intercept
        vals[21] = v_s3*(S4 - S3) + vals[20] + vals[19];// third intercept

    }


    else if (model_uid == 606)
    { 
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = params[3];
        double a_r = params[4];
        double a = params[5];
        double b = params[6];
        double m = params[7];
        double mt = params[8];
        double k1 = params[9];
        double k2 = params[10];
        double t_L = params[11];
        double NoFlow = params[12];
        double Td = params[13];
        double Beta = params[14];
        double lambda_1 = params[15];
        double lambda_2 = params[16];
        double v0 = params[17];

        vals[16] = 60.0*v0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[17] = v_r * (L_i / A_h) * 60; // [1/min] runoff speed.
    }

    else if (model_uid == 608)
    { 
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = params[3];
        double a_r = params[4];
        double a = params[5];
        double b = params[6]; 
        double c = params[7];
        double d = params[8];
        double k3 = params[9];
        double ki_fac = params[10];
        //double k2 = params[11];
        double t_L = params[11];
        double NoFlow = params[12];
        double Td = params[13];
        double Beta = params[14];
        double lambda_1 = params[15];
        double lambda_2 = params[16];
        double v0 = params[17];
        //double expo = params[18];

        vals[16] = 60.0*v0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[17] = v_r * (L_i / A_h) * 60; // [1/min] runoff speed.
    }

    else if (model_uid == 609)
    { 
         double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = params[3];
        double a_r = params[4];
        double a = params[5];
        double b = params[6]; 
        double c = params[7];
        double d = params[8];
        double k3 = params[9];
        double ki_fac = params[10];        
        double t_L = params[11];
        double NoFlow = params[12];
        double Td = params[13];
        double Beta = params[14];
        double lambda_1 = params[15];
        double lambda_2 = params[16];
        double v0 = params[17];
        //Pre computed parameters
        vals[16] = 60.0*v0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[17] = v_r * (L_i / A_h) * 60; // [1/min] runoff speed.
    }
    
    else if (model_uid == 610)
    { 
         double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = params[3];
        double a_r = params[4];
        double a = params[5];
        double b = params[6]; 
        double c = params[7];
        double d = params[8];
        double k3 = params[9];
        double ki_fac = params[10];        
        double t_L = params[11];
        double NoFlow = params[12];
        double Td = params[13];
        double Beta = params[14];
        double lambda_1 = params[15];
        double lambda_2 = params[16];
        double v0 = params[17];
        //double t_base = params[18];
        //Pre computed parameters
        vals[16] = 60.0*v0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[17] = v_r * (L_i / A_h) * 60; // [1/min] runoff speed.
    }

    else if (model_uid == 607)
    { 
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_r = params[3];
        double a_r = params[4];
        double a = params[5];
        double c = params[6];
        double m = params[7];
        double k1 = params[8];
        double k2 = params[9];
        double t_L = params[10];
        double Beta = params[11];
        double lambda_1 = params[12];
        double lambda_2 = params[13];
        double v0 = params[14];

        vals[13] = 60.0*v0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[14] = v_r * (L_i / A_h) * 60; // [1/min] runoff speed.
    }

    else if (model_uid == 654)
    {
        //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
        //The numbering is:     0   1   2    3     4   5   6   7 
        //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
        //The numbering is:        0      1        2     3   4     5        6   7  8 9    10     11
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_h = params[3];
        double k_i_factor = params[4];
        double k_3 = params[5];	//[1/min]
        double h_b = params[6];	//[m]
        double S_L = params[7];	//[m]
        double A = params[8];
        double B = params[9];
        double exponent = params[10];
        double lambda_1 = params[11];
        double lambda_2 = params[12]; 
        double v_0 = params[13];
               
       // vals[3] = v_h * L_i / A_h * 60.0; //Poner en el 3 el nuevo k_2   // [1/min] k_2
       // vals[4] = vals[3]* k_i_factor;                                  // [1/min] k_i
       // vals[12] = 60.0*v_0*pow(A_i, lambda_2)/((1.0 - lambda_1)*L_i);//[1/min]  invtau
       // vals[13] = A_h / 60.0;                                            // c_2
    }

    if (model_uid == 250)
    {
        //Order of parameters: A_i,L_i,A_h,h_r,invtau,k_2,k_I,c_1,c_2
        //The numbering is:	0   1   2   3    4     5   6   7   8
        //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,gamma,h_b,e_pot
        //The numbering is:        0      1        2     3   4     5         6    7	8
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double k_I_factor = global_params[5];

        vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
        vals[6] = vals[5] * k_I_factor;	//[1/min] k_I
        vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[8] = A_h / 60.0;	//  c_2
    }
    else if (model_uid == 252 || model_uid == 253)
    {
        //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
        //The numbering is:	0   1   2    3     4   5   6   7 
        //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
        //The numbering is:        0      1        2     3   4     5        6   7  8 9  10
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double k_i_factor = global_params[5];

        vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[4] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
        vals[5] = vals[4] * k_i_factor;	//[1/min] k_i
        vals[6] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[7] = A_h / 60.0;	//  c_2
    }
    else if (model_uid == 249 ||model_uid == 251 ||model_uid == 254 || model_uid == 256 || model_uid == 264)
    {
        //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
        //The numbering is:     0   1   2    3     4   5   6   7 
        //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
        //The numbering is:        0      1        2     3   4     5        6   7  8 9    10     11
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double k_i_factor = global_params[5];

        vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	 // [1/min]  invtau
        vals[4] = v_h * L_i / A_h * 60.0;                                // [1/min] k_2
        vals[5] = vals[4] * k_i_factor;                                  // [1/min] k_i
        vals[6] = (0.001 / 60.0);                                        // (mm/hr->m/min)  c_1
        vals[7] = A_h / 60.0;                                            // c_2
    }
    else if (model_uid == 255)
    {
        //Order of parameters: A_i,L_i,A_h,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent | invtau,k_2,k_i,c_1,c_2
        //The numbering is:	0   1   2   3   4      5       6   7  8 9   10        11    12  13  14  15
        //Order of global_params: v_0,lambda_1,lambda_2
        //The numbering is:        0      1        2
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_h = params[3];
        double k_i_factor = params[5];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];

        vals[11] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[12] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
        vals[13] = vals[12] * k_i_factor;	//[1/min] k_i
        vals[14] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[15] = A_h / 60.0;	//  c_2
    }
    else if (model_uid == 257)
    {
        //Order of parameters: A_i,L_i,A_h,horder,invtau,k_2,k_i,c_1,c_2
        //The numbering is:    0   1   2   3      4      5   6   7...8
        //Order of global_params: v_0_0,...,v_0_9,lambda_1_0,...,lambda_1_9,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A ,B, exponent,v_B,k_tl
        //The numbering is:       0         9     10,            19         20       21  22  23         24  25  26 27 28       29  30
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        int h_order = (int) params[3];
        assert(h_order < 10);

        double v_0 = global_params[h_order - 1];
        double lambda_1 = global_params[10 + h_order - 1];
        double lambda_2 = global_params[20];
        double v_h = global_params[21];
        double k_i_factor = global_params[23];

        vals[4] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[5] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
        vals[6] = vals[5] * k_i_factor;	//[1/min] k_i
        vals[7] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[8] = A_h / 60.0;	//  c_2
    }
    else if (model_uid == 258)
    {
        //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
        //The numbering is:	    0   1   2    3     4   5   6   7 
        //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I,h_b,S_L,v_B
        //The numbering is:        0      1        2     3   4   5   6   7   8
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double k_i_factor = global_params[5];

        vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
        vals[4] = v_h * L_i / A_h * 60.0;                               // [1/min] k_2
        vals[5] = vals[4] * k_i_factor;	                                // [1/min] k_i
        vals[6] = (0.001 / 60.0);                                       // (mm/hr->m/min)  c_1
        vals[7] = A_h / 60.0;                                           // c_2
    }
    else if (model_uid == 259)
    {
        //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
        //The numbering is:	    0   1   2    3     4   5   6   7 
        //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I,h_b,S_L,v_B,k_tl
        //The numbering is:        0      1        2     3   4   5   6   7   8   9
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        double v_h = global_params[3];
        double k_i_factor = global_params[5];

        vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
        vals[4] = v_h * L_i / A_h * 60.0;                               // [1/min] k_2
        vals[5] = vals[4] * k_i_factor;	                                // [1/min] k_i
        vals[6] = (0.001 / 60.0);                                       // (mm/hr->m/min)  c_1
        vals[7] = A_h / 60.0;                                           // c_2
    }
    else if (model_uid == 260)
    {
        //Order of parameters: A_i,L_i,A_h | invtau,c_1,c_2
        //The numbering is:	0   1   2  |    3    4   5 
        //Order of global_params: v_0,lambda_1,lambda_2,h_b,k_D,k_2,k_dry,k_i,T_L,N,phi
        //The numbering is:        0      1        2     3   4   5   6     7   8  9  10
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];

        vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[4] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[5] = A_h / 60.0;	//  c_2
    }
    else if (model_uid == 261)
    {
        //Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
        //The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
        //Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
        //The numbering is:        0      1        2    3  4   5 
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double S_h = params[3];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];

        vals[9] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[10] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[11] = A_h / 60.0;	//  c_2
        vals[12] = pow(S_h, 0.5)*L_i / A_h;	//c_3
    }
    else if (model_uid == 262)
    {
        //Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
        //The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
        //Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
        //The numbering is:        0      1        2    3  4   5
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double S_h = params[3];
        double eta = params[5];

        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];

        vals[10] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
        vals[11] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        vals[12] = A_h / 60.0;	//  c_2
        vals[13] = pow(S_h, 0.5)*L_i / (A_h*eta);	//k_2
    }
  else if (model_uid == 263) {
			double* vals = params;
			double A_i = params[0];
			double L_i = params[1];
			double A_h = params[2];
			double v_0=params[3];
			double lambda_1 = params[4];
			double lambda_2 = params[5];
			double v_h=params[6];
			double k_i_factor = params[7];
		    double k_3 = params[8];

		    double h_b = params[9];
		    double S_L = params[10];
		    double A = params[11];//global_params[8];
		    double B = params[12];//global_params[9];
		    double exponent = params[13];//global_params[10];
		    double v_B = params[14];//global_params[11];
		    double k_tl = params[15];//]global_params[12];


		}
	else if (model_uid == 300 || model_uid == 301) {
		/*
		 //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		 //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		 //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		 //The numbering is:        0      1        2     3   4   5
		 //Need to set entries 12-19 of params.
		 double* vals = params;
		 double K_T = 1.0;
		 double s_r = 1.0;
		 double rootS_h = pow(vals[7],.5);
		 double v_r = global_params[0];
		 double lambda_1 = global_params[1];
		 double lambda_2 = global_params[2];
		 double Q_r = global_params[3];
		 double A_r = global_params[4];
		 double RC = global_params[5];

		 vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);
		 vals[13] = vals[3] / s_r;
		 vals[14] = 2.0 / .6 * vals[0] * s_r / Q_r * 1.0 / vals[8] * rootS_h * pow(vals[3],2.0/3.0);
		 vals[15] = vals[6] * vals[0] * vals[3] / (3600.0 * Q_r);
		 vals[16] = 1e-3/(60*s_r) * RC;
		 vals[17] = 60 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * 1.0 / vals[8] * rootS_h;
		 vals[18] = K_T/60.0;
		 vals[19] = vals[6]/(60.0*s_r);

		 iparams[0] = link_i->location;
		 */
	} else if (model_uid == 315) {
		/*
		 //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		 //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		 //Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
		 //The numbering is:        0      1        2     3   4   5
		 //Need to set entries 12-19 of params.
		 double* vals = params;
		 double K_T = 1.0;
		 double s_r = 1.0;
		 double rootS_h = pow(vals[7],.5);
		 double L = params[0];
		 double A_h = params[1] * 1e6;	//Put into m^2
		 double eta = params[8];
		 double v_r = global_params[0];
		 double lambda_1 = global_params[1];
		 double lambda_2 = global_params[2];
		 double v_h = global_params[3];
		 double A_r = global_params[4];
		 double RC = global_params[5];


		 vals[12] = 60.0*v_r*pow(vals[2]/A_r,lambda_2)/((1.0-lambda_1)*vals[0]);	//invtau [1/min]
		 vals[13] = vals[3] / s_r; //epsilon
		 vals[14] = v_h*L;	//c_1 [m^2/s]
		 vals[15] = vals[6] * vals[0] * vals[3] / 3600.0; //c_2
		 vals[16] = (1e-3/60.0) * RC;	//c_3
		 vals[17] = 60.0*v_h*L/A_h;	//c_4 [1/min], A_h converted above
		 vals[18] = K_T/60.0;
		 vals[19] = vals[6]/(60.0*s_r);


		 iparams[0] = link_i->location;
		 */
	} else if (model_uid == 400) //tetis01 model
		{
		double* vals = params;
		double A_i = params[0]; //upstream area of the hillslope
		double L_i = params[1];	// channel lenght
		double A_h = params[2]; //area of the hillslope
		double v_0 = global_params[0]; //velocity river in channels [m/s]
		double lambda_1 = global_params[1]; //power discharge in routing function
		double lambda_2 = global_params[2]; //power of area in routing function
		double Hu = global_params[3]; //max available storage static storage [mm]
		double infiltration = global_params[4]; //infiltration rate [mm/hr]
		double percolation = global_params[5]; //percolation rate [mm/hr]
		double alfa2 = global_params[6]; //surface velocity [m/s]
		double alfa3 = global_params[7]; //linear reserv. coef gravitational storage [days]
		double alfa4 = global_params[8]; //linear reserv. coef aquifer storage [days]
        double melt_factor = global_params[9]; // melting factor in mm/hour/degree
        double temp_thres = global_params[10]; // in celsius degrees
		vals[3] = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L_i);//[1/min]  invtau params[3]
		vals[4] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
		vals[5] = A_h / 60.0;	//  c_2
    } else if (model_uid == 612) //tetis01 model
			{
		double* vals = params;
		double A_i = params[0]; //upstream area of the hillslope
		double L_i = params[1];	// channel lenght
		double A_h = params[2]; //area of the hillslope


		double v_0 = global_params[0]; //velocity river in channels [m/s]
		double lambda_1 = global_params[1]; //power discharge in routing function
		double lambda_2 = global_params[2]; //power of area in routing function
		double Hu = global_params[3]; //max available storage static storage [mm]
		double infiltration = global_params[4]; //infiltration rate [mm/hr]
		double percolation = global_params[5]; //percolation rate [mm/hr]
		double alfa2 = global_params[6]; //linear reservoir coef. surface storage [minutes]
		double alfa3 = global_params[7]; //linear reserv. coef gravitational storage [days]
		double alfa4 = global_params[8]; //linear reserv. coef aquifer storage [days]
		vals[3] = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L_i);//[1/min]  invtau params[3]
		vals[4] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
		vals[5] = A_h / 60.0;	//  c_2
        vals[6] = (L_i/A_h)*60; // c_3 [1/min] the residency time in the hillslope

	} else if (model_uid == 401 || model_uid == 402) //tetis02 & 03 model
			{
		double* vals = params;
		double A_i = params[0]; //upstream area of the hillslope
		double L_i = params[1];	// channel lenght
		double A_h = params[2]; //area of the hillslope
		double v_0 = global_params[0]; //velocity river in channels [m/s]
		double lambda_1 = global_params[1]; //power discharge in routing function
		double lambda_2 = global_params[2]; //power of area in routing function
		double Hu = global_params[3]; //max available storage static storage [mm]
		double infiltration = global_params[4]; //infiltration rate [mm/hr]
		double percolation = global_params[5]; //percolation rate [mm/hr]
		double alfa2 = global_params[6]; //surface velocity [m/s]
		double alfa3 = global_params[7]; //linear reserv. coef gravitational storage [days]
		double alfa4 = global_params[8]; //linear reserv. coef aquifer storage [days]
        double melt_factor = global_params[9]; // melting factor in mm/hour/degree
        double temp_thres = global_params[10]; // in celsius degrees
		vals[3] = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L_i);//[1/min]  invtau params[3]
		vals[4] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
		vals[5] = A_h / 60.0;	//  c_2

	} else if (model_uid == 2000) {
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		//Need to set entries 12-19 of params.
		double* vals = params;
		double K_T = 1.0;
		double s_r = 1.0;
		double rootS_h = pow(vals[7], .5);
		double v_r = global_params[0];
		double lambda_1 = global_params[1];
		double lambda_2 = global_params[2];
		double Q_r = global_params[3];
		double A_r = global_params[4];
		double RC = global_params[5];

		vals[12] = 60.0 * v_r * pow(vals[2] / A_r, lambda_2)
				/ ((1.0 - lambda_1) * vals[0]);
		vals[13] = vals[3] / s_r;
		vals[14] = 2.0 / .6 * vals[0] * s_r / Q_r * 1.0 / vals[8] * rootS_h
				* pow(vals[3], 2.0 / 3.0);
		vals[15] = vals[6] * vals[0] * vals[3] / (3600.0 * Q_r);
		vals[16] = 1e-3 / (60 * s_r) * RC;
		vals[17] = 60 * 1e-6 * 2.0 / .6 * vals[0] / vals[1] * 1.0 / vals[8]
				* rootS_h;
		vals[18] = K_T / 60.0;
		vals[19] = vals[6] / (60.0 * s_r);
	}
}

//Set the initial condition for the differential equations a link. This method will be called once for each link. The differential
//	variables from any .ini, .uini, or database are already set here. Precalculations have also been made previously. All
//	algebraic variables MUST be set here. Other initial conditions may be set here as well, depending on the model.
//VEC* global_params:	The vector of global parameters for the system.
//VEC* params:		The parameters for this link. DEM, dam parameters, and precalculations are available.
//unsigned int dam:	1 if a dam is present at this link, 0 if no dam is present.
//VEC* y_0:		The initial condition vector. Store the initial data here.
//unsigned int model_uid:	The index of the model.
//Returns the state of the solution (use 0 if state discontinuities are not a concern). !!!! Should the return value be an int? !!!!
int ReadInitData(
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    QVSData* qvs,
    unsigned short int dam,
    double *y_0, unsigned int dim,
    unsigned int model_uid,
    unsigned int diff_start, unsigned int no_init_start,
    void* user,
    void* external)
{
    unsigned int state;

    if (model_uid == 19)
    {
        //For model_uid 21, just set the state
        return LinearHillslope_Evap_Check(y_0, params, global_params, qvs, 0);
    }
    else if (model_uid == 21)
    {
        //For model_uid 21, only the storage (S or y_0[1]) has been set.
        //Order of parameters: A_i,L_i,A_h,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
        //The numbering is:	0   1   2  3  4    5	       6      7       8     9	  10	    11       12  13  14
        //Order of global_params: v_r,lambda_1,lambda_2,RC,S_0,v_h
        //The numbering is:        0      1        2     3  4   5
        double RC = global_params[3];
        double S_0 = global_params[4];
        double A_h = params[2];
        y_0[2] = RC * S_0 * A_h;
        y_0[3] = (1.0 - RC) * S_0 * A_h;

        state = dam_check(y_0, dim, global_params, num_global_params, params, num_params, qvs, dam, user);
        dam_q(y_0, dim, global_params, params, qvs, state, user, y_0);
        return state;
    }
    else if (model_uid == 22)
    {
        //For model_uid 22, only the storage (S or y_0[1]) has been set.
        //Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
        //The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
        //Order of global_params: lambda_1,lambda_2,S_0
        //The numbering is:         0        1       2
        double S_0 = global_params[2];
        double A_h = params[2];
        double RC = params[3];

        y_0[2] = RC * S_0 * A_h;
        y_0[3] = (1.0 - RC) * S_0 * A_h;

        state = dam_check2(y_0, dim, global_params, num_global_params, params, num_params, qvs, dam, user);
        dam_q2(y_0, dim, global_params, params, qvs, state, user, y_0);
        return state;
    }
    else if (model_uid == 23)
    {
        //For model_uid 23, only the storage (S or y_0[1]) has been set.
        //Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau,orifice_area,H_spill,H_max,S_max,alpha,orifice_diam,c_1,c_2,L_spill
        //The numbering is:	0   1   2  3   4   5   6  7   8          9	  10	  11   12     13      14        15  16   17
        //Order of global_params: lambda_1,lambda_2,S_0
        //The numbering is:         0        1       2
        double S_0 = global_params[2];
        double A_h = params[2];
        double RC = params[3];

        y_0[2] = RC * S_0 * A_h;
        y_0[3] = (1.0 - RC) * S_0 * A_h;

        state = dam_check3(y_0, dim, global_params, num_global_params, params, num_params, qvs, dam, user);
        dam_q3(y_0, dim, global_params, params, qvs, state, user, y_0);
        return state;
    }
    else if (model_uid == 40)
    {
        //For model_uid 40, only the storage (S or y_0[1]) has been set.
        //Order of parameters: A_i,L_i,A_h,RC,v_h,v_r,k2,k3,invtau
        //The numbering is:	0   1   2  3   4   5   6  7   8
        //Order of global_params: lambda_1,lambda_2,S_0,v_g
        //The numbering is:         0        1       2   3
        double S_0 = global_params[2];
        double A_h = params[2];
        double RC = params[3];

        y_0[2] = RC * S_0 * A_h;
        y_0[3] = (1.0 - RC) * S_0 * A_h;

        state = dam_check_qvs(y_0, dim, global_params, num_global_params, params, num_params, qvs, dam, user);
        dam_q_qvs(y_0, dim, global_params, params, qvs, state, user, y_0);

        return state;
    }
    else if (model_uid == 30)
    {
        //y_i has q (y_i[0]), s_p (y_i[1]), h_w (y_i[2]), theta (y_i[3])
        //Order of parameters: L_h,A_h,A_up,H_b,H_h,max_inf_rate,K_SAT,S_H,n_vh, b, c, d,K_Q,V_T,c_1,c_2,c_3,c_4,c_5,c_6,c_7
        //The numbering is:     0   1   2    3   4       5         6    7   8    9 10 11  12  13  14  15  16  17  18  19  20
        //Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,K_T,C_r,e_pot
        //The numbering is:        0      1        2     3   4   5   6   7
        /*
        double H_h = params[4];

        params[9] = 1.0;
        params[10] = 0.0;
        params[11] = 0.0;

        if(H_h < 1e-12)		//Flat surface
        {
        y_0[2] = 0.0;
        y_0[3] = 0.0;
        }
        */
        return 0;
    }
    else if ((model_uid == 191) || (model_uid == 192))
    {
        //For this model_uid, the extra states need to be set (3,4,5)
        y_0[3] = 0.0;
        y_0[4] = 0.0;
        y_0[5] = y_0[0];	//I'm not really sure what to use here...
    }
	else if (model_uid==193 || model_uid==194){return 0;}
    else if (model_uid == 195)
    {
        //For this model_uid, the extra states need to be set (3)
        y_0[3] = 0.0;
    }
    else if (model_uid == 196)
    {
        y_0[3] = 0.0;  // zero precip accumulation
        y_0[4] = 0.0;  // zero runoff accumulation
    } else if (model_uid == 200) {
		//For model_uid 200, only the discharge has been set. Need to set the storage.
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,RC,u_0
		//The numbering is:        0      1        2     3   4   5  6
		y_0[1] = params[0] / (global_params[6] + global_params[0]) * y_0[0];
		return 0;
	}else if (model_uid == 249) {
		//For this model_uid, the extra states need to be set (4,5,6)
		y_0[4] = y_0[0];
		y_0[5] = y_0[0];
	}
     else if (model_uid == 254) {
		//For this model_uid, the extra states need to be set (4,5,6)
		y_0[4] = 0.0;
		y_0[5] = 0.0;
		y_0[6] = y_0[0];
	} else if (model_uid == 255) {
		//Discharges are initially read into y_0[1] when no dam is present. So y_0[1] is copied to y_0[0],
		//then the corresponding storage is moved into y_0[1]. When a dam is present, y_0[1] will have the storage.
		//So the discharge can be calculated and stored into y_0[0].

		//Contains 2 layers in the channel: discharge, storage. Contains 3 layers on hillslope: ponded, top layer, soil.
		//Order of the states is:              0          1                                        2        3       4
		//Order of parameters: A_i,L_i,A_h,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent | invtau,k_2,k_i,c_1,c_2
		//The numbering is:	0   1   2   3   4      5       6   7  8 9   10        11    12  13  14  15
		//Order of global_params: v_0,lambda_1,lambda_2
		//The numbering is:        0      1        2

		if (dam) {
			unsigned int i;
			y_0[0] = y_0[1];
			for (i = 0; i < qvs->n_values - 1; i++)
				if (qvs->points[i][1] <= y_0[0]
						&& y_0[0] < qvs->points[i + 1][1])
					break;
			if (i == qvs->n_values - 1) {
				y_0[0] = qvs->points[i][1];
				y_0[1] = qvs->points[i][0];
			} else {
				double q2 = qvs->points[i + 1][1];
				double q1 = qvs->points[i][1];
				double S2 = qvs->points[i + 1][0];
				double S1 = qvs->points[i][0];
				y_0[1] = (S2 - S1) / (q2 - q1) * (y_0[0] - q1) + S1;
			}
			return i;
		} else {
			double lambda_1 = global_params[1];
			double tau_in_secs = 1.0 / params[11] * 60.0;
			y_0[0] = y_0[1];
			y_0[1] = tau_in_secs / (1.0 - lambda_1)
					* pow(y_0[0], 1.0 - lambda_1);
			return -1;
		}
	} else if (model_uid == 256 || model_uid == 264) {
		//For this model_uid, the extra states need to be set (4,5,6,7)
		y_0[4] = 0.0;
		y_0[5] = 0.0;
		y_0[6] = 0.0;
		y_0[7] = y_0[0];
	} else if (model_uid == 257) {
		//For this model_uid, the extra states need to be set (4,5,6,7)
		y_0[4] = 0.0;
		y_0[5] = 0.0;
		y_0[6] = 0.0;
		y_0[7] = y_0[0];
    }
    else if (model_uid == 200)
    {
        //For model_uid 200, only the discharge has been set. Need to set the storage.
        //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
        //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
        //Order of global_params: v_0,lambda_1,lambda_2,Q_r,A_r,RC,u_0
        //The numbering is:        0      1        2     3   4   5  6
        y_0[1] = params[0] / (global_params[6] + global_params[0]) * y_0[0];
        return 0;
    }
    else if (model_uid == 254)
    {
        //For this model_uid, the extra states need to be set (4,5,6)
        y_0[4] = 0.0;
        y_0[5] = 0.0;
        y_0[6] = y_0[0];
    }
    else if (model_uid == 255)
    {
        //Discharges are initially read into y_0[1] when no dam is present. So y_0[1] is copied to y_0[0],
        //then the corresponding storage is moved into y_0[1]. When a dam is present, y_0[1] will have the storage.
        //So the discharge can be calculated and stored into y_0[0].

        //Contains 2 layers in the channel: discharge, storage. Contains 3 layers on hillslope: ponded, top layer, soil.
        //Order of the states is:              0          1                                        2        3       4
        //Order of parameters: A_i,L_i,A_h,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent | invtau,k_2,k_i,c_1,c_2
        //The numbering is:	0   1   2   3   4      5       6   7  8 9   10        11    12  13  14  15
        //Order of global_params: v_0,lambda_1,lambda_2
        //The numbering is:        0      1        2

        if (dam)
        {
            unsigned int i;
            y_0[0] = y_0[1];
            for (i = 0; i<qvs->n_values - 1; i++)
                if (qvs->points[i][1] <= y_0[0] && y_0[0] < qvs->points[i + 1][1])	break;
            if (i == qvs->n_values - 1)
            {
                y_0[0] = qvs->points[i][1];
                y_0[1] = qvs->points[i][0];
            }
            else
            {
                double q2 = qvs->points[i + 1][1];
                double q1 = qvs->points[i][1];
                double S2 = qvs->points[i + 1][0];
                double S1 = qvs->points[i][0];
                y_0[1] = (S2 - S1) / (q2 - q1) * (y_0[0] - q1) + S1;
            }
            return i;
        }
        else
        {
            double lambda_1 = global_params[1];
            double tau_in_secs = 1.0 / params[11] * 60.0;
            y_0[0] = y_0[1];
            y_0[1] = tau_in_secs / (1.0 - lambda_1) * pow(y_0[0], 1.0 - lambda_1);
            return -1;
        }
    }
    else if (model_uid == 256)
    {
        //For this model_uid, the extra states need to be set (4,5,6,7)
        y_0[4] = 0.0;
        y_0[5] = 0.0;
        y_0[6] = 0.0;
        y_0[7] = y_0[0];
    }
    else if (model_uid == 257)
    {
        //For this model_uid, the extra states need to be set (4,5,6,7)
        y_0[4] = 0.0;
        y_0[5] = 0.0;
        y_0[6] = 0.0;
        y_0[7] = y_0[0];
    }
    else if ((model_uid == 258)||(model_uid == 259))
    {
        // For this model_uid, the extra states need to be set (4,5,6)
        y_0[4] = 0.0;
        y_0[5] = 0.0;
        y_0[6] = 0.0;
        y_0[7] = y_0[0];
    }
    else if (model_uid == 261)
    {
        //Discharges are initially read into y_0[1] when no dam is present. So y_0[1] is copied to y_0[0],
        //then the corresponding storage is moved into y_0[1]. When a dam is present, y_0[1] will have the storage.
        //So the discharge can be calculated and stored into y_0[0].

        //Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
        //The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
        //Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
        //The numbering is:        0      1        2    3  4   5 


        //For this model, the extra states need to be set (5,6,7)
        y_0[5] = 0.0;
        y_0[6] = 0.0;
        y_0[7] = y_0[0];

        if (dam)
        {
            int state = dam_check_qvs(y_0, dim, global_params, num_global_params, params, num_params, qvs, dam, user);
            dam_TopLayerNonlinearExpSoilvel(y_0, dim, global_params, params, qvs, state, user, y_0);
            return state;
        }
        else
        {
            double lambda_1 = global_params[1];
            double tau_in_secs = 1.0 / params[9] * 60.0;
            y_0[0] = y_0[1];
            //y_0[1] = pow(tau_in_secs*y_0[0],1.0-lambda_1);
            y_0[1] = tau_in_secs / (1.0 - lambda_1) * pow(y_0[0], 1.0 - lambda_1);
            return -1;
        }
    }
    else if (model_uid == 262)
    {
        //Discharges are initially read into y_0[1] when no dam is present. So y_0[1] is copied to y_0[0],
        //then the corresponding storage is moved into y_0[1]. When a dam is present, y_0[1] will have the storage.
        //So the discharge can be calculated and stored into y_0[0].

        //Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
        //The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
        //Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
        //The numbering is:        0      1        2    3  4   5

        //For this model, the extra states need to be set (5,6,7)
        y_0[5] = 0.0;
        y_0[6] = 0.0;
        y_0[7] = y_0[1];	//Note: See comment above. y_0[1] as the initial discharge.

        if (dam)
        {
            unsigned int i;
            y_0[0] = y_0[1];
            for (i = 0; i<qvs->n_values - 1; i++)
                if (qvs->points[i][1] <= y_0[0] && y_0[0] < qvs->points[i + 1][1])	break;
            if (i == qvs->n_values - 1)
            {
                y_0[0] = qvs->points[i][1];
                y_0[1] = qvs->points[i][0];
            }
            else
            {
                double q2 = qvs->points[i + 1][1];
                double q1 = qvs->points[i][1];
                double S2 = qvs->points[i + 1][0];
                double S1 = qvs->points[i][0];
                y_0[1] = (S2 - S1) / (q2 - q1) * (y_0[0] - q1) + S1;
            }
            return i;
        }
        else
        {
            double lambda_1 = global_params[1];
            double tau_in_secs = 1.0 / params[10] * 60.0;
            y_0[0] = y_0[1];
            y_0[1] = tau_in_secs / (1.0 - lambda_1) * pow(y_0[0], 1.0 - lambda_1);
            return -1;
        }
    }
    else if (model_uid == 300)
    {
        /*
        //For this model_uid, all initial conditions for variational equation must be set here.
        //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
        //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
        //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
        //The numbering is:        0      1        2     3   4   5
        unsigned int i;
        unsigned int offset = model_uid - 299;
        for(i=offset;i<y_0.dim;i++)	y_0[i] = 0.0;
        y_0[iparams[0] + offset] = 1.0;
        */
        return 0;
    }
    else if (model_uid == 301)
    {
        /*
        //For this model_uid, all initial conditions for variational equation must be set here.
        //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
        //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
        //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
        //The numbering is:        0      1        2     3   4   5
        unsigned int i;
        unsigned int offset = model_uid - 299;

        //New
        y_0[offset] = 1.0;
        y_0[offset + 1] = 1.0;
        y_0[offset + 2] = 0.0;
        for(i=offset+3;i<y_0.dim;i++)	y_0[i] = 0.0;
        */
        return 0;
    }
    else if (model_uid == 315)
    {
        //For this model_uid, all initial conditions for variational equation must be set here.
        //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
        //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
        //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
        //The numbering is:        0      1        2     3   4   5
        unsigned int i;
        unsigned int offset = 2;

        //New
        y_0[offset] = 1.0;
        y_0[offset + 1] = 1.0;
        y_0[offset + 2] = 0.0;
        for (i = offset + 3; i < dim; i++)
            y_0[i] = 0.0;

        return 0;
    } else if (model_uid == 261) {
		//Discharges are initially read into y_0[1] when no dam is present. So y_0[1] is copied to y_0[0],
		//then the corresponding storage is moved into y_0[1]. When a dam is present, y_0[1] will have the storage.
		//So the discharge can be calculated and stored into y_0[0].

		//Order of parameters: A_i,L_i,A_h,S_h,T_L,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,c_3
		//The numbering is:	0   1   2   3   4   5   6   7     8  |   9     10  11  12
		//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
		//The numbering is:        0      1        2    3  4   5

		//For this model, the extra states need to be set (5,6,7)
		y_0[5] = 0.0;
		y_0[6] = 0.0;
		y_0[7] = y_0[0];

		if (dam) {
			int state = dam_check_qvs(y_0, dim, global_params,
					num_global_params, params, num_params, qvs, dam, user);
			dam_TopLayerNonlinearExpSoilvel(y_0, dim, global_params, params,
					qvs, state, user, y_0);
			return state;
		} else {
			double lambda_1 = global_params[1];
			double tau_in_secs = 1.0 / params[9] * 60.0;
			y_0[0] = y_0[1];
			//y_0[1] = pow(tau_in_secs*y_0[0],1.0-lambda_1);
			y_0[1] = tau_in_secs / (1.0 - lambda_1)
					* pow(y_0[0], 1.0 - lambda_1);
			return -1;
		}
	} else if (model_uid == 262) {
		//Discharges are initially read into y_0[1] when no dam is present. So y_0[1] is copied to y_0[0],
		//then the corresponding storage is moved into y_0[1]. When a dam is present, y_0[1] will have the storage.
		//So the discharge can be calculated and stored into y_0[0].

		//Order of parameters: A_i,L_i,A_h,S_h,T_L,eta,h_b,k_D,k_dry,k_i | invtau,c_1,c_2,k_2
		//The numbering is:	0   1   2   3   4   5   6   7   8     9  |   10    11  12  13
		//Order of global_params: v_0,lambda_1,lambda_2,N,phi,v_B
		//The numbering is:        0      1        2    3  4   5

		//For this model, the extra states need to be set (5,6,7)
		y_0[5] = 0.0;
		y_0[6] = 0.0;
		y_0[7] = y_0[1];//Note: See comment above. y_0[1] as the initial discharge.

		if (dam) {
			unsigned int i;
			y_0[0] = y_0[1];
			for (i = 0; i < qvs->n_values - 1; i++)
				if (qvs->points[i][1] <= y_0[0]
						&& y_0[0] < qvs->points[i + 1][1])
					break;
			if (i == qvs->n_values - 1) {
				y_0[0] = qvs->points[i][1];
				y_0[1] = qvs->points[i][0];
			} else {
				double q2 = qvs->points[i + 1][1];
				double q1 = qvs->points[i][1];
				double S2 = qvs->points[i + 1][0];
				double S1 = qvs->points[i][0];
				y_0[1] = (S2 - S1) / (q2 - q1) * (y_0[0] - q1) + S1;
			}
			return i;
		} else {
			double lambda_1 = global_params[1];
			double tau_in_secs = 1.0 / params[10] * 60.0;
			y_0[0] = y_0[1];
			y_0[1] = tau_in_secs / (1.0 - lambda_1)
					* pow(y_0[0], 1.0 - lambda_1);
			return -1;
		}
	} else if (model_uid == 263) {
		//For this model_uid, the extra states need to be set (4,5,6,7)
		y_0[4] = 0.0;
		y_0[5] = 0.0;
		y_0[6] = 0.0;
		y_0[7] = y_0[0];
	}
	else if (model_uid == 300) {
		/*
		 //For this model_uid, all initial conditions for variational equation must be set here.
		 //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		 //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		 //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		 //The numbering is:        0      1        2     3   4   5
		 unsigned int i;
		 unsigned int offset = model_uid - 299;
		 for(i=offset;i<y_0.dim;i++)	y_0[i] = 0.0;
		 y_0[iparams[0] + offset] = 1.0;
		 */
		return 0;
	} else if (model_uid == 301) {
		/*
		 //For this model_uid, all initial conditions for variational equation must be set here.
		 //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		 //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		 //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		 //The numbering is:        0      1        2     3   4   5
		 unsigned int i;
		 unsigned int offset = model_uid - 299;

		 //New
		 y_0[offset] = 1.0;
		 y_0[offset + 1] = 1.0;
		 y_0[offset + 2] = 0.0;
		 for(i=offset+3;i<y_0.dim;i++)	y_0[i] = 0.0;
		 */
		return 0;
	} else if (model_uid == 315) {
		//For this model_uid, all initial conditions for variational equation must be set here.
		//Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
		//The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
		//Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
		//The numbering is:        0      1        2     3   4   5
		unsigned int i;
		unsigned int offset = 2;

		//New
		y_0[offset] = 1.0;
		y_0[offset + 1] = 1.0;
		y_0[offset + 2] = 0.0;
		for (i = offset + 3; i < dim; i++)
			y_0[i] = 0.0;

		return 0;
	} else if (model_uid == 400 || model_uid == 612)        //tetis01
		{

	    } 
    else if (model_uid == 401)        //tetis02
	    {

	    }
	else if (model_uid == 402)        //tetis03
			{
		//For this model_uid, the extra states need to be set (4,5,6,7)
		y_0[4] = 0.0;
		y_0[5] = 0.0;
		y_0[6] = 0.0;
		y_0[7] = y_0[0];
		//y_0[8] = 0.0;

//Discharges are initially read into y_0[1] when no dam is present. So y_0[1] is copied to y_0[0],
//then the corresponding storage is moved into y_0[1]. When a dam is present, y_0[1] will have the storage.
//So the discharge can be calculated and stored into y_0[0].
		if (dam) {
			unsigned int i;
			y_0[0] = y_0[8];
			for (i = 0; i < qvs->n_values - 1; i++)
				if (qvs->points[i][1] <= y_0[0]
						&& y_0[0] < qvs->points[i + 1][1])
					break;
			if (i == qvs->n_values - 1) {
				y_0[0] = qvs->points[i][1];
				y_0[8] = qvs->points[i][0];
			} else {
				double q2 = qvs->points[i + 1][1];
				double q1 = qvs->points[i][1];
				double S2 = qvs->points[i + 1][0];
				double S1 = qvs->points[i][0];
				y_0[8] = (S2 - S1) / (q2 - q1) * (y_0[0] - q1) + S1;
			}
			return i;
		} else {
			double lambda_1 = global_params[1];
			double tau_in_secs = 1.0 / params[10] * 60.0;
			y_0[0] = y_0[1];
			y_0[1] = tau_in_secs / (1.0 - lambda_1)
					* pow(y_0[0], 1.0 - lambda_1);
			return -1;
		}
	} else {
		//If not using algebraic variables, then everything is already set
		return 0;
	}

	return 0;
}

/*
 //If using data assimilation, sets the global errors and dim correctly
 void AssimError(unsigned int N,UnivVars* GlobalVars,ErrorData* GlobalErrors)
 {
 unsigned int i,j;
 unsigned int old_num_dense;

 //Remove any variaional indices from dense_indices
 for(i=0;i<globals->num_dense;i++)
 {
 if(globals->dense_indices[i] >= globals->problem_dim)
 {
 for(j=i+1;j<globals->num_dense;j++)	globals->dense_indices[j-1] = globals->dense_indices[j];
 (globals->num_dense)--;
 i--;
 }
 }

 old_num_dense = globals->num_dense;

 GlobalVars.dim = globals->problem_dim + N*globals->problem_dim*globals->problem_dim;
 globals->num_dense += N*globals->problem_dim*globals->problem_dim;
 GlobalErrors->abstol = realloc(GlobalErrors->abstol,GlobalVars.dim*sizeof(double));
 GlobalErrors->reltol = realloc(GlobalErrors->reltol,GlobalVars.dim*sizeof(double));
 GlobalErrors->abstol_dense = realloc(GlobalErrors->abstol_dense,GlobalVars.dim*sizeof(double));
 GlobalErrors->reltol_dense = realloc(GlobalErrors->reltol_dense,GlobalVars.dim*sizeof(double));
 GlobalErrors->abstol.dim = GlobalErrors->reltol.dim = GlobalErrors->reltol_dense.dim = GlobalErrors->reltol_dense.dim = GlobalVars.dim;

 //Setup error
 for(i=globals->problem_dim+1;i<GlobalVars.dim;i++)
 {
 GlobalErrors->abstol[i] = GlobalErrors->abstol[globals->problem_dim];
 GlobalErrors->reltol[i] = GlobalErrors->reltol[globals->problem_dim];
 GlobalErrors->abstol_dense[i] = GlobalErrors->abstol_dense[globals->problem_dim];
 GlobalErrors->reltol_dense[i] = GlobalErrors->reltol_dense[globals->problem_dim];
 }

 //Setup dense indices
 //Add in variational indices
 globals->dense_indices = realloc(globals->dense_indices,globals->num_dense * sizeof(unsigned int));
 for(i=old_num_dense;i<globals->num_dense;i++)
 {
 //globals->dense_indices[i] = i;
 globals->dense_indices[i] = (i-old_num_dense) + globals->problem_dim;
 }
 }
 */

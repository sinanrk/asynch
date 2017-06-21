Configuration File
==================

Configuration files are used to specify ALL inputs for the solvers. This includes river network topology files, model parameters, initial conditions, forcings, what information is printed to output files, etc...

Although each setting in a configuration file modifies values for the ASYNCH solvers, their exact use can be modified by the user altering the underlying source code. This can be done with calls to routines described in :ref:`C API`.

Global File Format
------------------

.. warning::
  As of version 1.4, the global file format is deprecated and may be removed in future releases. We encourage user to switch to the more robust and versatile :ref:`JSON File Format`.

Global files  (.gbl) are always ASCII files, assumed to be in UNIX format. Global files have a very rigid structure, unlike XML files, and information must be specified in a particular order. **The description of each input of the global files below are given in the order in which they are to be specified in the actual file**. The percent sign ``%`` is used for single line comments As described below, certain inputs are expected to be given within a single line Other than this restriction, white space is ignored Arguments surrounded by ``{ }`` below are mandatory, while those surrounded by the square brackets ``[ ]`` are dependent upon values specified by other parameters.

Overview
~~~~~~~~

Here is a typical global file taken from the examples folder:

::

  %Model UID
  190

  %Begin and end date time
  2017-01-01 00:00
  2017-01-02 00:00

  %Parameters to filenames
  0

  %Components to print
  3
  Time
  LinkID
  State0

  %Peakflow function
  Classic

  %Global parameters
  %6 v_r  lambda_1  lambda_2  RC    v_h  v_g
  6  0.33  0.20      -0.1     0.33  0.1  2.2917e-5

  %No. steps stored at each link and
  %Max no. steps transfered between procs
  %Discontinuity buffer size
  30 10 30

  %Topology (0 = .rvr, 1 = database)
  0 test.rvr

  %DEM Parameters (0 = .prm, 1 = database)
  0 test.prm

  %Initial state (0 = .ini, 1 = .uini, 2 = .rec, 3 = .dbc, 4 = .h5)
  1 test.uini

  %Forcings (0 = none, 1 = .str, 2 = binary, 3 = database, 4 = .ustr,
  %          5 = forecasting, 6 = .gz binary, 7 = recurring)
  2

  %Rain
  1 test.str

  %Evaporation
  7 evap.mon
  1398902400 1588291200

  %Dam (0 = no dam, 1 = .dam, 2 = .qvs)
  0

  %Reservoir ids (0 = no reservoirs, 1 = .rsv, 2 = .dbc file)
  0

  %Where to put write hydrographs
  %(0 = no output, 1 = .dat file, 2 = .csv file, 3 = database, 5 = .h5 packet, 6 = .h5 array)
  5 5.0 outputs.h5

  %Where to put peakflow data
  %(0 = no output, 1 = .pea file, 2 = database)
  1 test.pea

  %.sav files for hydrographs and peak file
  %(0 = save no data, 1 = .sav file, 2 = .dbc file, 3 = all links)
  1 test.sav %Hydrographs
  3 %Peakflows

  %Snapshot information (0 = none, 1 = .rec, 2 = database, 3 = .h5, 4 = recurrent .h5)
  4 60 test.h5

  %Filename for scratch work
  tmp

  %Numerical solver settings follow

  %facmin, facmax, fac
  .1 10.0 .9

  %Solver flag (0 = data below, 1 = .rkd)
  0
  %Numerical solver index (0-3 explicit, 4 implicit)
  2
  %Error tolerances (abs, rel, abs dense, rel dense)
  1e-3 1e-3 1e-3
  1e-6 1e-6 1e-6
  1e-3 1e-3 1e-3
  1e-6 1e-6 1e-6

  # %End of file
  -------------------------------

In the following sections, we will go into details about the meaning and options for each entry.

Model Type
~~~~~~~~~~

Format:

::

  {model id}

This value specifies the id for the model to be used. This is a non-negative integer value which corresponds to a particular system of ordinary-differential equations (or possibly DAEs). Examples of built-in models is given in  :ref:`Built-in Models`. If you are using the API to use  custom models, this model id is ignored, see :ref:`Custom Models`.

Simulation period
~~~~~~~~~~~~~~~~~

Format:

::

  {begin datetime}
  {end datetime}

The begin and end datetimes are given in ``YYYY-MM-DD HH:MM`` format using the UTC timezone or in `Unix Time`_ format.

Parameters on Filenames
~~~~~~~~~~~~~~~~~~~~~~~

Format:

::

  {parameter on output filename flag}

This is a boolean value (``0`` or ``1``) that indicates whether all output filenames should include the uniform in space and time parameters ``0`` indicates no, ``1`` indicates yes. This feature can be useful for keeping track of output files from multiple simulations.

Solver Outputs
~~~~~~~~~~~~~~

Format:

::

  {number of outputs}
  [output1]
  [output2]
  [...]

This set of input parameters specifies the names of all outputs from the solvers. Several built in outputs exist, and the user is able to construct his own outputs. Built in outputs are given in :ref:`Built-In Output Time Series`. Output names are case sensitive. The first required value is the number of outputs (>= 0), followed by the names of each output, on separate lines.

Peakflow Statistics Function Name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Format:

::

  {function name}

This sets the function for outputting peakflow information. The built in peakflow function "Classic" is one option, and the user is free to construct his own. A function name must be specified here, even if peakflow information is not requested for any links.

Global Parameters
~~~~~~~~~~~~~~~~~

Format:

::

  {number of parameters} [parameter 1] [parameter 2] ... [parameter n]

This is where model parameters which are constant in space and time are specified. The first value is a nonnegative integer specifying the number of global parameters to follow. Every model requires a certain number of global parameters. If the number given in the global file is less than expected for a particular model, an error occurs. If the number is greater than expected, a warning is given. These "extra" parameters are available to the model for use. This can sometimes be useful for quick tests, but should be avoided normally.

The parameter meanings depend upon the model used. The units of these parameters is also model dependent.

Buffer Sizes
~~~~~~~~~~~~

Format:

::

  {steps stored at each link} {max number of steps transferred} {discontinuity buffer size}

These nonnegative integer values allow the user to specify sizes of internal buffers. In general, as these numbers are increased, the solvers run faster, but more memory is required. A good starting point that works in most cases is the set ``30 10 30``. Typically, if these values need to be less than 10 to run the solvers, a deeper issue with memory constraints should be addressed.

Topology
~~~~~~~~

Format:

::

  {topology flag} [output link id] {.rvr filename or .dbc filename}

This is where connectivity of the river network is specified. This can be done in one of two ways If the topology flag is ``0``, a river topology file (.rvr) is used. If the topology flag is ``1``, then topology is downloaded from the database specified with the database file (.dbc). The database connection allows for one additional feature: a subbasin can be specified If the output link id is taken to be 0, all link ids found in the database are used. Otherwise, the link with link id specified and all upstream links are used. Pulling subbasins from a topology file is not currently supported.

Link Parameters
~~~~~~~~~~~~~~~

Format:

::

  {parameter flag} {.prm filename or .dbc filename}

This specifies where parameters which vary by link and not time, are specified If the parameter flag is ``0``, the parameters are given in a parameter (.prm) file. If the flag is ``1``, then the parameters are downloaded from the database specified by the database connection file (.dbc). The number, order, meaning, and units of these parameters varies from model to model.

Initial States
~~~~~~~~~~~~~~

Format:

::

  {initial state flag} {.ini, .uini, .rec, .dbc or .h5 filename} [unix time]

This section specifies the initial state of the model. The values for the initial state flag can be ``0``, ``1``, ``2``, ``3`` or ``4`` corresponding, respectively, to a ini, uini, rec, dbc, h5 file. The unix time argument is used for database connections only. This value is available in the query of the database connection file and can be used for selecting values from tables.

Forcings
~~~~~~~~

Format:

::

  {number of forcings}
  [forcing1 flag] [forcing1 information]
  [forcing2 flag] [forcing2 information]
  [...]

Information about time dependent forcings is specified here. Each model has an expected number of forcings. If the number of forcings specified here is less than expected, an error is thrown. If the number of forcings is greater than expected, a warning is given. This warning allows for tests to be performed and implemented quickly. In general, this feature should be avoided.

Forcing information varies considerably based upon the corresponding forcing flag. Several forcing types require unix times to determine what forcing data to use. If a model requires multiple forcings with unix times, the times do not need to be consistent, i.e., one forcing could start on July 1st 2014 at midnight, while another forcing starts at April 5th 2008.

No Forcing
^^^^^^^^^^

Format:

::

  0

A forcing flag of ``0`` specifies no forcing input. This is the same as a forcing value of `0.0` for all links and all time.

Storm File
^^^^^^^^^^

Format:

::

  1 {.str filename}

A forcing flag of ``1`` indicates the forcing is specified by a .str file. The filename and path of a valid storm (.str) file is required.

Binary Files
^^^^^^^^^^^^

Format:

::

  2 {binary file identifier}
  {chunk size} {time resolution} {beginning file index} {ending file index}

A forcing flag of ``2`` indicates the forcing is specified by a collection of binary forcing files. The identifier can be adorned with a path to the binary files. The chunk size is a positive integer that indicates the number of binary files kept in memory at once. The time resolution indicates the amount of time between successively indexed binary files. This value is a floating point number with units equal to those of the time variable of the model used The beginning and ending file indices indicate the range of the binary files to be used. The indices are integer valued. The simulation will begin using the binary file with index given by the beginning file index. If the total simulation time would require binary files with index greater than the ending file index, the forcing values are taken to be 0.0 for all such binary files.

Forcings from Databases
^^^^^^^^^^^^^^^^^^^^^^^

Format:

::

  3 {.dbc filename}
  {chunk size} {time resolution} {beginning unix time} {ending unix time}

A forcing flag of ``3`` indicates the forcing data will be pulled from a PostgreSQL database. The database connection filename can include a path. The chunk size is a positive integer representing the number of forcing values pulled from the database at once from each link. A chunk size of 10 tends to work well. A larger chunk size requires more memory and larger datasets returned from the database, but a small number of queries. The time resolution is a floating point number with units in minutes. This represents the time resolution of the data in the accessed table. The integrity of the database table is not thoroughly checked by the solvers.

The simulation will begin using the data from the database with unix time given by the beginning unix time. If the total simulation time would require data from the database with unix time greater than the ending unix time, the forcing values are taken to be 0.0 for times greater than the ending unix time.

Uniform Forcings
^^^^^^^^^^^^^^^^

Format:

::

  4 {.ustr filename}

A forcing flag of ``4`` indicates a forcing that is uniform in space. The forcings are given by a uniform storm file (.ustr).

GZipped Binary
^^^^^^^^^^^^^^

Format:

::

  6 {gzipped binary file identifier}
  {chunk size} {time resolution} {beginning file index} {ending file index}

A forcing flag of ``6`` indicates the forcing is specified by a collection of binary forcing files that have been gzipped (compressed as .gz files). All parameters for this flag are identical to that of using binary files with forcing flag ``3``.

Monthly Forcings
^^^^^^^^^^^^^^^^

Format:

::

  7 { mon filename}
  {beginning unix time} {ending unix time}

A forcing flag of ``7`` indicates a uniform in space forcing that recurs monthly. When the end of the calendar year is reached, the monthly forcing file (.mon) is read again from the beginning The beginning unix time is used to determine the month the simulation begins (for this forcing). If the total simulation time takes the simulation past the ending unix time, the forcing is assumed to be ``0.0`` for all locations and times beyond the ending unix time

Grid Cell
^^^^^^^^^

Format:

::

  8 {index filename}
  {chunk size} {beginning file index} {ending file index}

A forcing flag of ``8`` indicates the forcing is specified by a collection of grid cell forcing files. The index filename can be adorned with a path to the index file. The chunk size is a positive integer that indicates the number of grid cell files kept in memory at once. The beginning and ending file indices indicate the range of the grid cell files to be used. The indices are integer valued.

The simulation will begin using the grid cell file with index given by the beginning file index. If the total simulation time would require grid cell files with index greater than the ending file index, the forcing values are taken to be ``0.0`` for all such grid cell files. In addition, if a grid cell file is missing, all values at each cell are assumed to be ``0.0``.

Dams
~~~~

Format:

::

  {dam flag} [.dam or .qvs filename]

This section specifies whether dams will be used A dam flag of ``0`` means no dams are used. A flag of ``1`` indicates a dam file ( dam) will be used, and a flag value of ``2`` indicates a discharge vs storage file ( qvs) will be used. Some models do not support dams. For these models, the dam flag must be set to ``0`` or an error occurs.

State Forcing Feeds
~~~~~~~~~~~~~~~~~~~

Format:

::

  {reservoir flag} [.rsv or .dbc filename] [forcing index]

This section specifies whether a provided forcing (see :ref:`Forcings`) is to be used as a forcing of the states of differential or algebraic equations at some links. A reservoir flag of ``0`` indicates no forcing will by applied to system states. A flag of ``1`` indicates state forcings will be applied to all link ids in the specified .rsv file. A reservoir flag of ``2`` indicates state forcing will be applied to all link ids pulled from the database the given .dbc file. If the reservoir flag is not ``0``, then the index of the forcing must be specified.

Time Series Location
~~~~~~~~~~~~~~~~~~~~

Format:

::

  {time series flag} [time resolution] [.dat or .csv or .h5 or .dbc filename] [table name]

This section specifies where the final output time series will be saved. A time series flag value of ``0`` indicates no time series data will be produced. Any flag with value greater than ``0`` requires a time resolution for the data. This value has units equal to the units of total simulation time (typically minutes). A value of ``-1`` uses a resolution which varies from link to link based upon the expression:

.. math::

  \begin{align}
   \left(0.1 \cdot \frac{A}{1 \ km^2} \right)^{\frac{1}{2}} \ min
  \end{align}

where :math:`A` is the upstream of the link, measured in km2.

A time series flag of ``1`` indicates the results of the simulation will be saved as a .dat file. The filename complete with a path must be specified. If a file with the name and path given already exists, it is overwritten.
A time series flag of ``2`` indicates the results will be stored as a .csv file.
A time series flag of ``3`` indicates the results will be uploaded into the database described by the given .dbc file. In this case, a table name accessible by the queries in the .dbc file must be specified.
A time series flag of ``5`` indicates the results will be stored as a .h5 HDF5 file with a packet layout compatible with PyTable.
A time series flag of ``6`` indicates the results will be stored as a .h5 HDF5 file with an 3D array layout. Time, link id and output indexes are given as additional 1D "dimension" arrays. Selected outputs in :ref:`Solver Outputs` must have the same type (ASYNCH_FLOAT).

This section is independent of the section for Link IDs to Save described below (see :ref:`Global Parameters`) For example, if link ids are specified in the Link IDs to Save section and the time series flag in the Time Series Locations set to ``0``, no output is generated. Similarly, if *the time series id flag* is set to ``0`` in the Link IDs to Save section and the time series flag is set to ``1``, a .dat file with ``0`` time series is produced.

.. note::

  The time resolution is entirely independent of the time step used by the numerical integrators. Reducing this value does NOT produce more accurate results. To improve accuracy, reduce the error tolerances described in :ref:`Numerical Error Tolerances`. There is no built-in way to produce results at every time step, as this is a very easy way to crash a compute node or file system.

Peakflow Data Location
~~~~~~~~~~~~~~~~~~~~~~

Format:

::

  {peakflow flag} [.pea / .dbc filename] [table name]

This section specifies where the final peakflow output will be saved. A peakflow flag of ``0`` indicates no peakflow data is produced. A peakflow flag of ``1`` indicates the peakflow results of the simulation will be saved as a .pea file. The filename complete with a path from the binary file must be specified. A peakflow flag of ``2`` indicates the results will be uploaded into the database described by the given .dbc file. In this case, a table name accessible by the queries in the dbc file must be specified.

This section is independent of the section for Link IDs to Save described below (see :ref:`Link IDs to Save`). For example, if link ids are specified in the Link IDs to Save section and the peakflow flag in the peakflow data location is set to ``0``, no output is generated. Similarly, if the peakflow id flag is set to ``0`` in the Link IDs to Save section and the peakflow flag is set to ``1``, a .pea file with ``0`` peakflows is produced.

Link IDs to Save
~~~~~~~~~~~~~~~~

Format:

::

  {time series id flag} [.sav / .dbc filename]
  {peakflow id flag} [.sav / .dbc filename]

This section provides the list of link ids in which data is produced. The first line is for the time series outputs, while the second is for the peakflow outputs. The time series ID flag and the peakflow ID flag take the same list of possible values. A flag of ``0`` indicates no link IDs for which to produce data. A flag of ``1`` indicates the list of link IDs is provided by the corresponding save file (.sav). A flag of ``2`` indicates the list of link IDs is provided by the database specified in the given database connection file (.dbc). A flag of ``3`` indicates that all links will have data outputted.

.. warning::

  A time series ID flag of ``3`` can easily wreak havoc on a file system for simulations with a large number of links. At the very least, extremely large output files and database tables will occur. Be very careful with this! Typically, using a flag value of ``3`` for peakflow link ids, or for the time series ID flag for a very small basin (< 500 links) will not create any problems.

This section is independent of the sections for Time Series Location and peakflow data location above (see :ref:`Time Series Location` and :ref:`Peakflow Data Location`). For example, if link ids are specified in the Link IDs to Save section and the time series flag in the Time Series Location set to ``0``, no output is generated. Similarly, if the time series id flag is set to ``0`` in the Link IDs to Save section and the time series flag is set to ``1``, a .dat file with zero time series is produced.

Snapshot Information
~~~~~~~~~~~~~~~~~~~~

Format:

::

  {snapshot flag} [time step of periodical snapshots] [.rec / .dbc / .h5 filename]

This section specifies where snapshot information is produced. A snapshot is a record of *every state* at *every link* in the network. Snapshots can be produced at the end of simulations or periodically. This is useful for beginning a new simulation where an old one ended. A snapshot flag of ``0`` indicates no snapshot is produced. A snapshot flag of ``1`` indicates the snapshot will be produced as a recovery (.rec) file with path and filename specified. A snapshot flag of ``2`` indicates the snapshot will be uploaded to the database specified by the database connectivity (.dbc) file.

A snapshot flag of ``3`` indicates the snapshot will be produced as a HDF5 (.h5) file with path and filename specified. A snapshot flag of ``4`` generates periodical snapshots in which case an addition parameter gives the interval between two snapshots and the second parameter is the output basename. For example:

::

  %Snapshot information (0 = none, 1 = .rec, 2 = .dbc, 3 = .h5, 4 = periodical .h5)
  4 60 filename.h5

generates

::

  filename_1480000000.h5
  filename_1480003600.h5
  filename_1480007200.h5
  ...

Scratch Work Location
~~~~~~~~~~~~~~~~~~~~~

Format:

::

  {filename}

This section specifies the location of temporary files. These files are used to store intermediate calculations. The filename can include a path name. If the file already exists, the contents are overwritten. If a simulation is aborted, these files may not be removed. Otherwise, they are deleted at the end of the simulation.

Error Control Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

Format:

::

  {facmin} {facmax} {fac}

This section specifies parameters related to the error control strategy of the numerical integrators. The value facmin represents the largest allowed decrease in the stepsize of the integrators as a percent of the current step Similarly, facmax represents the largest allowed increase. The value fac represents the safety factor of the integrators. Any accepted stepsize is multiplied by this value Good values of facmin, facmax, and fac to use are ``0`` 1, 10 0, and ``0`` 9, respectively

Numerical Error Tolerances
~~~~~~~~~~~~~~~~~~~~~~~~~~

Format:

::

  {solver flag} [ rkd filename]
  [rk solver index]
  [absolute error tolerance 1] [absolute error tolerance 2]
  [relative error tolerance 1] [relative error tolerance 2]
  [dense absolute error tolerance 1] [dense absolute error tolerance 2]
  [dense relative error tolerance 1] [dense relative error tolerance 2]

This section specifies error tolerances for the numerical integrators. A solver flag of ``0`` indicates the same tolerances will be used for all links. A solver flag of ``1`` indicates the tolerance info will be specified in the given RK data (.rkd) file. If solver flag is ``0``, than an rk solver index must be specified. A list of Runge-Kutta methods is given in :ref:`Built-In Runge-Kutta Methods`. Each error tolerance must have a value for each state of the system. The order of the tolerances must match the order of the states in the state vectors. The absolute and relative error tolerances are those typically used for RK methods. The dense tolerances are for the numerical solution produced between time steps. A numerical solution is rejected if either the error tolerances or dense error tolerances for any state is believed to be violated.


JSON File Format
----------------

JavaScript Object Notation or JSON is an open-standard format that uses human-readable text to transmit data objects consisting of attribute–value pairs and array data types (or any other serializable value).

Overview
~~~~~~~~

Here is a typical JSON file for the :ref:`Constant Runoff Hydrological Model` taken from the examples folder:

.. code-block:: javascript

  {
    "model": 190,                 //Model UID
    "begin": "2017-01-01 00:00",  //Can be either YYYY-MM-DD HH:MM
    "end": 1483336800,            //or Unix Time
    "outputs": {
      "functions": [
        "Time",
        "LinkID",
        "State0"
      ],
      "timeseries": {
        "filename": "outputs.h5",
        "locations": "hydrographs.sav",
        "interval": 60
      },
      "peaks": {
        "filename": "peaks.pea",
        "locations": "peaks.sav",
        "function": "Classic"
      },
      "postfix": false,
    },
    "global_params": [ 0.33, 0.20, -0.1, 0.33, 0.1, 2.2917e-5 ],
    "topology": "test.rvr",
    "local_params": "test.prm",
    "initial_state": "test.uini",
    "forcings": [
      "test.str",
      {
        "filename": "evap.mon",
        "begin": 1398902400,
        "end": 1588291200
      }
    ],
    "solver": {
      "method": 2,
        "error_ctl" {
        "facmin": 0.1,
        "facmax": 10.0,
        "fac": 0.9
      },
      "tolerances": [
        [ 1e-3, 1e-6, 1e-3, 1e-6 ],
        [ 1e-3, 1e-6, 1e-3, 1e-6 ],
        [ 1e-3, 1e-6, 1e-3, 1e-6 ]
      ],
      "buffers": {
        "num_step": 30,
        "num_transfer": 10,
        "num_discont": 30
      }
    }
  }

.. note::  Comments are not allowed by the `JSON standard <http://www.json.org>`_.

Top-level properties
~~~~~~~~~~~~~~~~~~~~

Properties are required unless otherwise designated.

+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| Property      | Type                                          | Required? | Description                                                                                                      |
+===============+===============================================+===========+==================================================================================================================+
| model         | Integer                                       | Yes       | Model Unique Identifier                                                                                          |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| begin         | Integer or String                             | Yes       | The begin datetime is given in YYYY-MM-DD HH:MM format using the UTC timezone or in `Unix Time`_ format.         |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| end           | Integer or String                             | Yes       | The end datetime is given in YYYY-MM-DD HH:MM format using the UTC timezone or in `Unix Time`_ format.           |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| outputs       | :ref:`Object <Output properties>`             | No        | The outputs configuration. If omitted, no output is generated.                                                   |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| global_params | Array                                         | Yes       | The global parameters.                                                                                           |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| topology      | String                                        | Yes       | Path of the the river network topology file (.rvr or .dbc).                                                      |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| local_params  | String                                        | Yes       | Path of the local parameters file (.prm or .dbc).                                                                |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| initial_state | String                                        | Yes       | Path of the initial state of the system file (.ini, .uini, .rec, .dbc or .h5).                                   |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| forcings      | :ref:`Array <Forcing properties>`             | No        | The forcings configuration. If omitted, the model runs without forcing.                                          |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| snaphosts     | String or :ref:`Object <Snapshot properties>` | No        | Path of the snapshots file (.rec, .dbc or .h5) or snapshots configuration. If omitted, no snapshot is generated. |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+
| solver        | :ref:`Object <Solver properties>`             | No        | The solver configuration. If omitted, the :ref:`Default <Solver properties>` configuration is used.              |
+---------------+-----------------------------------------------+-----------+------------------------------------------------------------------------------------------------------------------+


Output properties
~~~~~~~~~~~~~~~~~

+------------+--------------------------------------+-----------+-----------+---------------------------------------------------------------------------+
| Property   | Type                                 | Required? | Default   | Description                                                               |
+============+======================================+===========+===========+===========================================================================+
| functions  | Array                                | Yes       | None      | The names of all outputs from the solvers.                                |
+------------+--------------------------------------+-----------+-----------+---------------------------------------------------------------------------+
| timeseries | :ref:`Object <Timeserie properties>` | No        | None      | The timeserie configuration. If omitted no timeserie output is generated. |
+------------+--------------------------------------+-----------+-----------+---------------------------------------------------------------------------+
| peaks      | :ref:`Object <Peak properties>`      | No        | None      | The peak configuration. If omitted no peak outut is generated.            |
+------------+--------------------------------------+-----------+-----------+---------------------------------------------------------------------------+
| postfix    | Boolean                              | No        | ``false`` | Postfix output filenames with the global parameters.                      |
+------------+--------------------------------------+-----------+-----------+---------------------------------------------------------------------------+

Timeserie properties
^^^^^^^^^^^^^^^^^^^^

+-----------+---------+-----------+-----------------------------------------------------------------------------------+
| Property  | Type    | Required? | Description                                                                       |
+===========+=========+===========+===================================================================================+
| filename  | String  | Yes       | The path of the timeserie output file.                                            |
+-----------+---------+-----------+-----------------------------------------------------------------------------------+
| locations | String  | No        | The path of the list of link ids to save file. If omitted, every links are saved. |
+-----------+---------+-----------+-----------------------------------------------------------------------------------+
| interval  | Integer | Yes       | The interval between two consecutive outputs.                                     |
+-----------+---------+-----------+-----------------------------------------------------------------------------------+

Peak properties
^^^^^^^^^^^^^^^

+-----------+---------+-----------+-------------+-----------------------------------------------------------------------------------+
| Property  | Type    | Required? | Default     | Description                                                                       |
+===========+=========+===========+=============+===================================================================================+
| filename  | String  | Yes       | None        | The path of the peak output file.                                                 |
+-----------+---------+-----------+-------------+-----------------------------------------------------------------------------------+
| locations | String  | Yes       | None        | The path of the list of link ids to save file. If omitted, every links are saved. |
+-----------+---------+-----------+-------------+-----------------------------------------------------------------------------------+
| function  | String  | No        | ``Classic`` |  The peakflow function name (see `Built-In Peakflow Functions`_)                  |
+-----------+---------+-----------+-------------+-----------------------------------------------------------------------------------+

Forcing properties
~~~~~~~~~~~~~~~~~~

+------------+---------+-----------+---------+------------------------------------------------------+
| Property   | Type    | Required? | Default | Description                                          |
+============+=========+===========+=========+======================================================+
| filename   | String  | Yes       | None    |                                                      |
+------------+---------+-----------+---------+------------------------------------------------------+
| chunk_size | Integer | Yes       | ``10``  | The number of forcing values kept in memory at once. |
+------------+---------+-----------+---------+------------------------------------------------------+
| time_step  | Number  | Yes       | None    | The time resolution.                                 |
+------------+---------+-----------+---------+------------------------------------------------------+

Snapshot properties
~~~~~~~~~~~~~~~~~~~

+----------+--------+-------------------------------------------------+
| Property | Type   | Description                                     |
+==========+========+=================================================+
| filename | String | The path of the snapshot file (.h5).            |
+----------+--------+-------------------------------------------------+
| interval | Number | The interval between two consecutive snapshots. |
+----------+--------+-------------------------------------------------+

Solver properties
~~~~~~~~~~~~~~~~~

The default configuration is:

.. code-block:: javascript

  {
    "method": 2,
    "error_ctl": {
      "facmin": 0.1,
      "facmax": 10.0,
      "fac": 0.9
    },
    "tolerances": [
      [ 1e-3, 1e-6, 1e-3, 1e-6 ],
      ...
    ],
    "buffers": {
      "num_step": 30,
      "num_transfer": 10,
      "num_discont": 30
    }
  }

All properties are optionals and have default values.

+----------------+------------------------------------------+----------------------------------------------+--------------------------------------------------------------------------------+
| Property       | Type                                     | Default                                      | Description                                                                    |
+================+==========================================+==============================================+================================================================================+
| method         | Integer                                  | ``2``                                        | A Runge-Kutta method (see :ref:`Built-In Runge-Kutta Methods`).                |
+----------------+------------------------------------------+----------------------------------------------+--------------------------------------------------------------------------------+
| error_ctl      | :ref:`Object <Error Control properties>` | :ref:`Default <Error Control properties>`    | Parameters related to the error control strategy of the numerical integrators. |
+----------------+------------------------------------------+----------------------------------------------+--------------------------------------------------------------------------------+
| tolerances     | :ref:`Array <Error Tolerance array>`     | :ref:`Default <Error Tolerance array>`       | Error tolerances for the numerical integrators for every state variable.       |
+----------------+------------------------------------------+----------------------------------------------+--------------------------------------------------------------------------------+
| buffers        | :ref:`Object <Buffer properties>`        | :ref:`Default <Buffer properties>`           | The sizes of internal buffers.                                                 |
+----------------+------------------------------------------+----------------------------------------------+--------------------------------------------------------------------------------+
| scratch_folder | String                                   | ``/tmp/``                                    | The location of temporary files.                                               |
+----------------+------------------------------------------+----------------------------------------------+--------------------------------------------------------------------------------+

Error Control properties
^^^^^^^^^^^^^^^^^^^^^^^^

+----------+--------+----------+---------------------------------------------------------------------------------------------------+
| Property | Type   | Default  | Description                                                                                       |
+==========+========+==========+===================================================================================================+
| facmin   | Number | ``0.1``  | The largest allowed decrease in the stepsize of the integrators as a percent of the current step. |
+----------+--------+----------+---------------------------------------------------------------------------------------------------+
| facmax   | Number | ``10.0`` | The largest allowed increase in the stepsize of the integrators as a percent of the current step. |
+----------+--------+----------+---------------------------------------------------------------------------------------------------+
| fac      | Number | ``0.9``  | The safety factor of the integrators.                                                             |
+----------+--------+----------+---------------------------------------------------------------------------------------------------+

Error Tolerance array
^^^^^^^^^^^^^^^^^^^^^

+-------+--------+----------+-------------------------------------------------------+
| Index | Type   | Default  | Description                                           |
+=======+========+==========+=======================================================+
| [0]   | Number | ``1e-3`` | The absolute error tolerance used for RK methods.     |
+-------+--------+----------+-------------------------------------------------------+
| [1]   | Number | ``1e-6`` | The relative error tolerance used for RK methods.     |
+-------+--------+----------+-------------------------------------------------------+
| [2]   | Number | ``1e-3`` | The absolute tolerances used by dense output methods. |
+-------+--------+----------+-------------------------------------------------------+
| [3]   | Number | ``1e-6`` | The relative tolerances used by dense output methods. |
+-------+--------+----------+-------------------------------------------------------+

Buffer properties
^^^^^^^^^^^^^^^^^

+--------------+---------+---------+----------------------------------------------------+
| Property     | Type    | Default | Description                                        |
+==============+=========+=========+====================================================+
| num_step     | Integer | ``30``  | Number of dense outputs steps stored at each link. |
+--------------+---------+---------+----------------------------------------------------+
| num_transfer | Integer | ``10``  | Number of dense outputs transferred at once.       |
+--------------+---------+---------+----------------------------------------------------+
| num_discont  | Integer | ``30``  | Number of discontinuity buffers at each link.      |
+--------------+---------+---------+----------------------------------------------------+

.. _`Unix Time`: https://en.wikipedia.org/wiki/Unix_time

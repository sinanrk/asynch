Input/Output Formats
====================

Numerous inputs are required to construct and solve the model equations. This includes link connectivity, link parameters, numerical error tolerances, precipitation forcings, a mathematical model, etc. Most inputs are specified on a link-by-link basis, and can be given as either a file or a database table.

In addition, ASYNCH must know information about data outputs. This includes where to store the outputs, what the outputs are, at which links outputs will be given, etc.

The format for each input and output is described below. For database inputs/outputs, ASYNCH can access PostgreSQL through the libpq libraries. Obviously, the user must have read or write permissions to a PostgreSQL database to utilize database features.

Database Connection Files
-------------------------

Database connection files are ASCII text files with a .dbc extension which specify how to connect to a database, and the queries to pull/push data from/to the database. Although the format of database connection files is the same, the specific requirements of the queries varies with how the file is used. For instance, queries for pulling link connectivity information is very different from queries for uploading peakflows to a database table.

Format:

::

  dbname={db} host={host} user={user} password={pass}
  {number of queries}
  [query 1]
  [query 2]
  ...

The first line of every database connection file specifies the information needed to make a connection. A user must have permission to read or write from the database at the given host; otherwise, queries sent to the database will fail. The number of queries will vary depending upon how the database connection file is used. The appropriate number of queries and what they should return is documented in the remainder of :ref:`Input/Output Formats`. The number of queries may be zero.

Any queries listed **MUST** be ended with a semicolon (;). For some queries, further information outside the database connection file may be available, depending upon how the query is to be used. This additional information is explained in the appropriate section below for input formats. Such information includes link ids and unix times. To denote in a query where this information should be placed, use the symbol ``"%u"`` for integers and ``"%s"`` for names.

Link Connectivity Input
-----------------------

Link connectivity information is used to specify how links in a network are connected. Topology can be provided through either a river network file (.rvr) file or through a database table. When a river network file is used, every link in the file is used (i.e. no subnetworks) then pulling connectivity data from a database, a subset of the network can be used.

Regardless of the format, all link ids must be given in the same order in the link connectivity, link parameter, and initial state inputs.

Rvr Files
~~~~~~~~~

River network files are ASCII text files with the following format:

::

  {number of links}
  {link id 1}
  {number of parents} [parent id 1] [parent id 2]
  {link id 2}
  {number of parents} [parent id 1] [parent id 2]

White space can be used freely throughout the file. The layout in the above specification is purely optional; the order of the information is what is important. The file begins with the total number of links in the file. Then each link id is specified, followed by the number of parents for the link and each of their ids. A link id can appear in a list of parent link ids at most once. If a link does not have parents, it must still appear in this file with a ``0`` for the number of parents.

Topology Database Queries
~~~~~~~~~~~~~~~~~~~~~~~~~

If the connectivity is pulled from a database, a corresponding database connection file is used This file requires three queries:

1. Query to pull all link ids from a table

  -  Inputs: none
  -  Returned tuples: (link id)

2. Query to pull all link id, parent link id pairs

  -  Inputs: none
  -  Returned tuples: (link id, parent link id)

3. Query to pull all link id, parent link id pairs upstream from a given outlet link id

  -  Inputs: outlet link id
  -  Returned tuples: (link id, parent link id)

The last two queries must return a tuple for each link id for each parent link. So a link with two parents should appear twice in the returned tuples, once for each parent link. The returned tuples must be grouped by the link id so all parent information appears consecutively.

Link Parameter Input
--------------------

Link parameter input specifies the parameters for the model that vary link to link This information can be provided in a parameter file ( prm) or through a database table The number of parameters for each link, their meaning, and their order depends upon the model used In particular, the value of disk params determines the number of parameters expected at each link See :ref:`SetParamSizes`.

Regardless of the format, all link ids must be given in the same order in the link connectivity, link parameter, and initial state inputs.

Prm File
~~~~~~~~

A parameter file is an ASCII text file with the following format:

::

  {number of links}
  {link id 1} {parameter 1} {parameter 2} {parameter 3}
  {link id 2} {parameter 1} {parameter 2} {parameter 3}

White space can be used freely throughout the file. The layout in the above specification is purely optional; the order of the information is what is important. The file begins with the total number of links. Then each link id is specified, followed by the parameters for that link.

Param Database Queries
~~~~~~~~~~~~~~~~~~~~~~

If the parameters are pulled from a database, a corresponding database connection file is used This file requires two queries:

1.  Query to pull all parameters

  -  Inputs: none
  -  Returned tuples: (link id, parameter 1, parameter 2, ...)

2.  Query to pull all parameters above an outlet

  -  Inputs: outlet link id
  -  Returned tuples: (link id, parameter 1, parameter 2, ...)

Initial Values Input
--------------------

The link initial values input specifies the initial values for the states of the differential and algebraic model equations This information can be provided in several different formats: an initial value file (.ini), a uniform initial value file (.uini), a recovery file (.rec), and through a database table.

Ini Files
~~~~~~~~~

An initial value file is an ASCII text file that lists the initial values for each link. The format is:

::

  {model type}
  {number of links}
  {initial time}
  {link id 1}
  {initial value 1} {initial value 2}
  {link id 2}
  {initial value 1} {initial value 2}

The model type is the number of the model to be used. This determines how many initial values are expected for the model. Initial states must be provided only for those states determined by differential equations, and only for those which require an initial condition. These are the states with index between ``diff_start`` and ``no_ini_start`` in the state vectors See :ref:`SetParamSizes`.

Uini Files
~~~~~~~~~~

A uniform initial value file is similar to an initial value file, but the initial values, when required, are the same at every link The format is given by:

::

  {model type}
  {initial time}
  {initial value 1} {initial value 2}

The model type is the number of the model to be used. This determines how many initial values are expected for the model. Initial values must be provided only for those states determined by differential equations, and only for those which require an initial condition. These are the states with index between ``diff_start`` and ``no_ini_start`` in the state vectors. See :ref:`SetParamSizes`. Notice that unlike an initial value file, no link ids are given, and only one set of initial values are given.

Rec Files
~~~~~~~~~

A recovery file is an ASCII text file that lists the initial values for each link. The format is:

::

  {model type}
  {number of links}
  {initial time}
  {link id 1}
  {initial value 1} {initial value 2}
  {link id 2}
  {initial value 1} {initial value 2}

The format is identical to that of an initial value file, with one important exception The initial value of EVERY state must be provided at each link. For models with ``diff_start`` set to 0 and ``no_ini_start`` set to dim, a recovery file is identical to an initial value file. See :ref:`SetParamSizes`.

.. warning:: For the initial values of algebraic equations, no checks on the input data are performed to ensure the solution is consistent.

Ini Database Queries
~~~~~~~~~~~~~~~~~~~~

If the initial values are pulled from a database, a corresponding database connection file is used. This file requires one query:

1. Query to pull all initial states for every link:

  -  Inputs: integer value
  -  Returned tuples: (link id, initial value 1, initial value 2, )

The query allows for one input to be used to obtain the needed information. This value could be, for example, an outlet link id or a unix time. Similar to recovery files, initial values must be provided for every link.

Ini HDF5 Files
~~~~~~~~~~~~~~

H5 Ini files are H5 that contains a single resizable Packet Tables or PyTable `snapshot`. Two HDF5 tools can be used to get the structure of the snapshots, ``l5hs`` and ``h5dump``:

.. code-block:: sh

  >h5ls -v test_1483228800.h5
  snapshot                 Dataset {11/Inf}
      Location:  1:1024
      Links:     1
      Chunks:    {512} 14336 bytes
      Storage:   308 logical bytes, 80 allocated bytes, 385.00% utilization
      Filter-0:  deflate-1 OPT {5}
      Type:      struct {
                     "link_id"          +0    native unsigned int
                     "state_0"          +4    native double
                     "state_1"          +12   native double
                     "state_2"          +20   native double
                 } 28 bytes

.. code-block:: sh

  >h5dump -H test_1483228800.h5
  HDF5 "test_1483228800.h5" {
  GROUP "/" {
     ATTRIBUTE "model" {
        DATATYPE  H5T_STD_U16LE
        DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
     }
     ATTRIBUTE "unix_time" {
        DATATYPE  H5T_STD_U32LE
        DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
     }
     ATTRIBUTE "version" {
        DATATYPE  H5T_STRING {
           STRSIZE 4;
           STRPAD H5T_STR_NULLTERM;
           CSET H5T_CSET_ASCII;
           CTYPE H5T_C_S1;
        }
        DATASPACE  SCALAR
     }
     DATASET "snapshot" {
        DATATYPE  H5T_COMPOUND {
           H5T_STD_U32LE "link_id";
           H5T_IEEE_F64LE "state_0";
           H5T_IEEE_F64LE "state_1";
           H5T_IEEE_F64LE "state_2";
        }
        DATASPACE  SIMPLE { ( 11 ) / ( H5S_UNLIMITED ) }
     }
  }
  }

Three global attributes are available :

+------------+--------------------------------------------------+
| Name       | Description                                      |
+============+==================================================+
| version    | The version of ASYNCH used to generate this file |
+------------+--------------------------------------------------+
| model      | The model id used to generate this file          |
+------------+--------------------------------------------------+
| issue_time | The unix time at the beginning of the time serie |
+------------+--------------------------------------------------+

Forcing Inputs
--------------

Numerous and diverse formats are implemented for feeding forcing inputs into a model. These formats vary considerably, and can have different impacts on performance.

Storm Files
~~~~~~~~~~~

Storm files (.str) provide an ASCII text format for setting a forcing at each link. The format of these files is:

::

  {number of links}
  {link id 1} {number of changes}
  {time 1} {value 1}
  {time 2} {value 2}
  {time 3} {value 3}
  {link id 2} {number of changes}
  {time 1} {value 1}
  {time 2} {value 2}
  {time 3} {value 3}

The format requires a time series to be provided for every link. The number of values can vary from link to link, and the time steps do not need to be uniformly spaced or even the same for each link. The first time at each link must be the same however, and correspond to the beginning of the simulation (typically 0). The forcings are assumed to be constant between time steps. After the final time step, the forcing value is held at the last value for the remainder of the simulation. The data provided by a storm file is entirely read into memory at the beginning of a run. As such, this format may not be suitable for large or long simulations.

Uniform Storm Files
~~~~~~~~~~~~~~~~~~~

Uniform storm files (.ustr) provide an ASCII text format for setting a forcing uniformly at each link. The format of these files is:

::

  {number of changes}
  {time 1} {value 1}
  {time 2} {value 2}
  {time 3} {value 3}

The format is similar to that of a storm file, but only one time series is specified, and is applied at every link. The time steps do not need to be uniformly spaced. The first time must correspond to the beginning of the simulation (typically 0) The forcing is assumed to be constant between time steps. After the fnal time step, the forcing value is held at the last value for the remainder of the simulation. The data provided by a uniform storm file is entirely read into memory at the beginning of a run. As such, this format may not be suitable for extremely long simulations.

Binary Storm Files
~~~~~~~~~~~~~~~~~~

Binary storm files (no extension) may also be used. Instead of a single file, these are a collection of files providing forcing values at different times. The format of these files is:

::

  {link id 1 value}
  {link id 2 value}
  {link id 3 value}

Each file is simply a list of forcing values for each link. Because link ids are not present, the values are assumed to be in the same order as the link ids from the link connectivity input. Further, a value must be specified for every link. The filename of each binary file should end with an integer value. The index of each file should be numbered consecutively. The time step size between files, as well as the range of files to use, are specified in the global file (see 6 1 10). If the simulation goes beyond the last file, all further forcing values are assumed to be ``0`` for all links. The values in binary storm files are stored as single precision floating point numbers, with endianness different from the native machine. These files are read into memory chunk by chunk. This allows for a ceiling on the memory used, independent of the number of files.

Gzipped Binary Storm Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

Gzipped binary storm files (.gz) are also supported. The format of these files is identical to that of the binary storm files, as they are simply gzipped binary storm files. The files are uncompressed before use, making this format slower than the regular binary storm files.

Grid Cell Files
~~~~~~~~~~~~~~~

Grid cell files group link ids together, and specify forcing values for each group (referred to as a cell). Although similar to binary files, this format requires two accompanying text files: an index file and a lookup file.

The index file specifies meta data for the grid cell data. The format for this ASCII file is:

::

  {time resolution (mins)}
  {conversion factor}
  {number of cells}
  {grid cell data file name prefix}
  {lookup filename}

The time resolution specifies the amount of time between grid cell files. The resolution is typically given in minutes. The conversion factor is a floating point number. Each forcing value in the grid cell files is multiplied by this factor. The number of cells is an integer value. Each grid cell filename has a prefix, followed by a number. The prefix is specified in the index file. The prefix may include a path. If the path is relative (i e , does not begin with a '/'), the path is taken relative to the location of the index file. Lastly, the index file includes the filename for the lookup file. A path is allowed, but is taken relative to the location of the index file, unless given as an absolute path.

The lookup file specifies which link ids are located in each cell. The format for this text file is:

::

  {link id 1} {cell id 1}
  {link id 2} {cell id 2}

The cell ids are indexed starting from 0, and the cell index cannot be larger than the number of cells specified in the accompanying index file.

The grid cell files are binary files. Each gives forcing values for each cell at a moment of time. If a cell is omitted from a file, then the forcing value is assumed to be 0. The filename for each grid cell file begins with the same prefix, followed by an integer. This integer indicated the order in which the grid cell files should be used. Although the number on these files should be consecutive, a missing file indicates all cells take a forcing value of 0. The format of the binary grid cell files is:

::

  {1} {cell id 1} {forcing value 1}
  {cell id 2} {forcing value 2}

The grid cell files are binary, so the spacing above is purely for readability. Each file begins with the integer value 1, stored as a 32-bit integer. This is used for checking the file format. Each cell id is given as a 32-bit integer and each forcing value is given as a 16-bit integer. Before the forcing values are loaded into memory, they are multiplied by the conversion factor given in the index file. Again, every cell does not need to be given in every grid cell file; only when the forcing value is nonzero does a value need to be given.

Monthly Recurring Forcing
~~~~~~~~~~~~~~~~~~~~~~~~~

Monthly recurring forcing files (.mon) allow forcings to be set monthly. These files are given as ASCII text files in the format:

::

  {value for January}
  {value for February}
  {value for March}
  ...
  {value for December}

A value must be given for each month. Over each simulated month, the forcing value is held constant, and is uniform over all links.

Database Forcing
~~~~~~~~~~~~~~~~

If the forcing data is pulled from a database, a corresponding database connection file is used. This file requires three queries:

1. Query to pull rainfall data for all link ids present in the database table

  -  Inputs: lower unix time, upper unix time
  - Returned tuples: (unix time, forcing value, link id)

2. Query to pull rainfall data for all link ids upstream from an outlet link

  - Inputs: outlet link id, lower unix time, upper unix time
  - Returned tuples: (unix time, forcing value, link id)

3. Query to pull a single forcing intensity from the database table

  -  Inputs: none
  -  Returned tuple: (unix time)

The first and second queries are similar, except that the second pulls only a subset of links from the database table. Forcing values are read in blocks from the table, with forcing values coming between the lower (inclusive) and upper (exclusive) unix times available to the queries. If a value for a link is not pulled from the database, the value is assumed to be 0.

The last query is used to find an actual valid timestamp in the database table. It needs to return only one unix time.

Dam Parameters Input
--------------------

Two formats currently exist for setting parameters at links with dams: dam parameter files (.dam) and discharge vs storage files (.qvs).

Dam Files
~~~~~~~~~

The format of dam parameter files is similar to that of parameter files:

::

  {number of links with dams}
  {link id 1}
  {parameter 1} {parameter 2} {parameter 3} ...
  {link id 2}
  {parameter 1} {parameter 2} {parameter 3} ...
  ...

The number of parameters needed for each link is model dependent and determined by the value dam params size. See :ref:`SetParamSizes`. For dam parameter files, only the links with dams must be listed here. Only links with id appearing in this file will have dams.

QVS Files
~~~~~~~~~

Discharge vs storage files take a series of discharge values and a corresponding series of storage values to decide the relationship between two states. The format of these files is similar to storm files (see :ref:`Forcing Inputs`):

::

  {number of links with dams}
  {link id 1} {number of points for link 1}
  {storage value 1} {discharge value 1}
  {storage value 2} {discharge value 2}
  ...
  {link id 2} {number of points for link 2}
  {storage value 1} {discharge value 1}
  {storage value 2} {discharge value 2}
  ...

The number of points at each link can vary For dam parameter files, only links with dams are listed here. Only links with id appearing in this file will have dams. Internally, discharges and storages with values between two points are interpolated. This interpolation process is model dependent.

Time Series Output
------------------

Three formats are supported for outputting time series calculations: data files (.dat), comma-separated values (.csv), and a database table. The particular time series calculated is set in the global file (see :ref:`Time Series Location`). The structure of each format is considerably different.

Data Files
~~~~~~~~~~

Data files are in ASCII text format. These files are designed to be generic and flexible so as to be easily read by whatever data analysis program the user prefers. Data files are created with the format:

::

  {number of links}
  {number of output values}
  {link id 1} {number of points for link id 1}
  {value 1 for series 1} {value 1 for series 2} {value 1 for series 3} ...
  {value 2 for series 1} {value 2 for series 2} {value 2 for series 3} ...
  ...
  {link id 2} {number of points for link id 2}
  {value 1 for series 1} {value 1 for series 2} {value 1 for series 3}
  {value 2 for series 1} {value 2 for series 2} {value 2 for series 3}
  ...

The series for the links appear in a column The number of points can vary from link to link, depending upon the user's selection in the global file The number of output values determines how many values appear in each line of the time series.

CSV Files
~~~~~~~~~

A CSV file is a typical format to make data easy to read in spreadsheet software. The structure of CSV files is:

::

  {link id 1},, ... , {link id 2},
  Output 1, Output 2,  , Output 1, Output 2,
  {value 1,1,1}, {value 1,2,1},  , {value 1,1,2}, {value 1,2,2},
  {value 2,1,1}, {value 2,2,1},  , {value 2,1,2}, {value 2,2,2},
  {value 3,1,1}, {value 3,2,1},  , {value 3,1,2}, {value 3,2,2},

The series for the links appear in a row. Under link id 1, each requested series appears, followed by the series for link id 2, and so on.

Out Database Queries
~~~~~~~~~~~~~~~~~~~~

A database connection file can be used to upload results into a database table This file requires only one query:

1. Query to create a table for uploading data

  - Inputs: table name
  - Returned tuples: none

The query should create the table where the series information is to be stored ASYNCH does NOT remove any existing data from the table, or check if the table exists already.

.. _out-hdf5-files:

Out HDF5 Files
~~~~~~~~~~~~~~

H5 outputs files are the prefered output format as it is both compact and efficient. There are also easy to read with third party software, see :ref:`Reading the HDF5 outputs with Python` for example.

H5 Output files are H5 that contains a single resizable Packet Tables or PyTable `outputs`. Two HDF5 tools can be used to get the strucutre of the outputs, ``l5hs`` and ``h5dump``:

.. code-block:: sh

  > l5hs -v outputs.h5
  Opened "outputs.h5" with sec2 driver.
  outputs                  Dataset {578/Inf}
      Location:  1:1024
      Links:     1
      Chunks:    {512} 8192 bytes
      Storage:   9248 logical bytes, 3455 allocated bytes, 267.67% utilization
      Filter-0:  deflate-1 OPT {5}
      Type:      struct {
                     "Time"             +0    native double
                     "LinkID"           +8    native int
                     "State0"           +12   native float
                 } 16 bytes

.. code-block:: sh

  >h5dump -H  outputs.h5
  HDF5 "outputs.h5" {
  GROUP "/" {
     ATTRIBUTE "issue_time" {
        DATATYPE  H5T_STD_U32LE
        DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
     }
     ATTRIBUTE "model" {
        DATATYPE  H5T_STD_U16LE
        DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
     }
     ATTRIBUTE "version" {
        DATATYPE  H5T_STRING {
           STRSIZE 4;
           STRPAD H5T_STR_NULLTERM;
           CSET H5T_CSET_ASCII;
           CTYPE H5T_C_S1;
        }
        DATASPACE  SCALAR
     }
     DATASET "outputs" {
        DATATYPE  H5T_COMPOUND {
           H5T_IEEE_F64LE "Time";
           H5T_STD_I32LE "LinkID";
           H5T_IEEE_F32LE "State0";
        }
        DATASPACE  SIMPLE { ( 578 ) / ( H5S_UNLIMITED ) }
     }
  }
  }

Three global attributes are available :

+------------+--------------------------------------------------+
| Name       | Description                                      |
+============+==================================================+
| version    | The version of ASYNCH used to generate this file |
+------------+--------------------------------------------------+
| model      | The model id used to generate this file          |
+------------+--------------------------------------------------+
| issue_time | The unix time at the beginning of the time serie |
+------------+--------------------------------------------------+


Peakflow Output
---------------

Peakflow outputs can be created in two formats: peakflow files (.pea) and database tables.

PEA Files
~~~~~~~~~

Peakflow files created with the "Classic" peakflow function take the structure:

::

  {number of link}
  {model type}
  {link id 1} {peakflow value} {time to peak} {area}
  {link id 2} {peakflow value} {time to peak} {area}
  {link id 3} {peakflow value} {time to peak} {area}

The time to peak is measured since the beginning of the simulation. The peakflow value for each link is the maximum value achieved over the simulation for the state with index ``0`` in the state vector. The area given is the parameter from the link parameters with index area idx. See :ref:`SetParamSizes`.

Peakflow Database Queries
~~~~~~~~~~~~~~~~~~~~~~~~~

Peakflow output may be written to a database table if a database connection file is specified. One query is required, and one additional query is optional:

1. Query to create a table for uploading data

  - Inputs: table name
  - Returned tuples: none

2. Query to delete contents from a table

  - Inputs: table name
  - Returned tuples: none

The first query should create the table where the peakflow information is to be stored ASYNCH does NOT remove any existing data from the table, or check if the table exists already. The second query is optional, and will be used to delete any existing contents from the table before uploading data. The particular values uploaded to the database are determined through the peakflow function defined in :ref:`Peakflow Statistics Function Name`.

Link IDs for Time Series and Peakflows
--------------------------------------

Link ids must be specified for time series output and peakflow outputs. This can be done in one of two formats: save files (.sav) and database tables. Each of these formats is effectively just a list of link ids.

SAV Files
~~~~~~~~~

The structure of save files is:

::

  {link id 1}
  {link id 2}
  {link id 3}
  ...

If a link id is specified in a save file, but is not present in the network, a warning will be issued, and the link id is ignored.

Save Database Queries
~~~~~~~~~~~~~~~~~~~~~

For pulling links from a database table, only one query is required:

1.  Query to pull link ids

  -  Inputs: none
  -  Returned tuples: (link id)

Snapshot Output
---------------

Snapshot outputs can take multiple formats: files and database tables. The format for recovery and hdf5 files is covered in :ref:`Initial Values Input` as an input.

For using a database table, a database connection file is specified. The database connection file has three optional queries:

1. Query to create a database table before the first upload

  -  Inputs: table name
  -  Returned tuples: none

2.  Query to delete a table before the first upload

  -  Inputs: table name
  -  Returned tuples: none

3.  Query to truncate a table before every upload

  -  Inputs: table name
  -  Returned tuples: none

In practice, snapshots are often applied periodically as a way to create check points for the program. The third query allows the user to limit the number of snapshots in a table to one.

Runge-Kutta Data Input
----------------------

Runge-Kutta Data files (.rkd) allow information about the underlying numerical methods to be specified link by link. These files are ASCII. The structure is given by:

::

  {number of links}
  {number of states}
  {link id 1}
  [absolute error tolerance 1] [absolute error tolerance 2]
  [relative error tolerance 1] [relative error tolerance 2]
  [dense absolute error tolerance 1] [dense absolute error tolerance 2]
  [dense relative error tolerance 1] [dense relative error tolerance 2]
  {RK Index for link id 1}
  {link id 2}
  [absolute error tolerance 1] [absolute error tolerance 2]
  [relative error tolerance 1] [relative error tolerance 2]
  [dense absolute error tolerance 1] [dense absolute error tolerance 2]
  [dense relative error tolerance 1] [dense relative error tolerance 2]
  {RK Index for link id 2}
  ...

An error tolerance is specified for every state at every link. The order of the links must match with the order given by the topology input, and number of states must agree with what the model expects.

Temporary Files
---------------

In general, sufficient memory is not available to hold a history of states while further computations take place. Thus, temporary files are created by ASYNCH to store time series outputs. These files are aggregated (typically after the simulation has completed) into the final time series output files (see :ref:`Time Series Output`).

Most users will not need to concern themselves with the underlying structure of these files. However, some routines exist for manipulating how these files are used, and an understanding of the temporary file structure is necessary for using those routines.

Temporary files are always in binary format. Every MPI process has its own file, which contains time series outputs for links assigned to the process. The format of the data looks like:

::

  {link id 1} {expected number of steps}
  {step 0, value 0} {step 0, value 1}
  {step 1, value 0} {step 1, value 1}
  {step 2, value 0} {step 2, value 1}
  {link id 2} {expected number of steps}
  {step 0, value 0} {step 0, value 1}
  {step 1, value 0} {step 1, value 1}
  {step 2, value 0} {step 2, value 1}
  ...

Because these files are a binary format, the presentation above is simply for readability. No new lines or spaces are present in the actual files. Only link ids for which time series data has been requested will appear in the files. Before any data is written, dummy values are placed in the files when they are created to insure the files are large enough. The number of these dummy files for each link is given by the expected number of steps value. This number is determined based upon the values of ``maxtime`` and the time resolution of the output time series when the temporary files are created.

.. warning::

   Modifcations to these values after creation of the temporary files could create a situation where the files are not large enough to contain every point in the time series. This generally happens when ``maxtime`` is increased. Therefore, when the temporary files are created, the largest expected value of ``maxtime`` should be set. If the temporary files are not large enough, a warning will be displayed during calculations, and further outputs will be lost.

While calculations are performed, processes will write output data to these files, replacing whatever dummy values are present. Modifying the behavior of these files is generally not needed, but can be performed through various routines. See :ref:`C API`. The link ids and expected number of steps in temporary files are represented by unsigned 32 bit integers. The data types and sizes for the time series data will vary depending upon the desired outputs. If a built-in output time series is used for the states, these will appear in the temporary files as 64 bit floating point numbers.

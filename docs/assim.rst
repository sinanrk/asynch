Data Assimilation
=================

In this section, we discuss the current implementation of the Data Assimilation and the usage of the ``assim`` command line interface.

Introduction
------------

The idea is to combine our hydrological model with stream flow observations.

Introduction to Data Assimilation and 3D-Var from my slides.

Variational equations
~~~~~~~~~~~~~~~~~~~~~

system state sensitivities

The 3D-Var algorithm is a simplification of the full variational data assimilation scheme, making the assumption that model field is stationary within a wide time-range called assimilation window. In that case, the observation within the window although depending on time, are considered as observation of the analysis time.


Least square fitting
~~~~~~~~~~~~~~~~~~~~




Building
--------

Data Assimilation is available only if asynch is built with the ``PETSc`` library. Refer to the :ref:`Installation` for more information. Make sure that ``./configure`` returns:

::

  checking for PETSC... yes

Usage
-----

``assim`` requires an  additional  configuration ``.das`` file on the command line, for exemple:

.. code-block:: sh

  assim turkey_river.gbl turkey_river.das

Configuration
-------------

Overview
~~~~~~~~

Here is a typical ``.das`` file taken from the examples folder:

::

  %Model variant
  254_q

  %Observation input (discharge)
  %Time resolution of observations
  assim51.dbc 15.0

  %Step to use (assimilation window)
  %12 % use 3 hours
  %24 % use 6 hours
  48 % use 12 hours
  %96 % use 24 hours

  %Max least squares iterations
  5

  # %End of file
  -------------------------------

Model variant
~~~~~~~~~~~~~

Format:

::

  {model id}

This string value specifies which assimilation model is used and which state variable initial conditions are optimized.

======= =============== ===
Id      Model           State variable
======= =============== ===
254     Top Layer Model Every state variable
254_q   Top Layer Model Discharge
254_qsp Top Layer Model Discharge, pond storage
254_qst Top Layer Model Discharge, top layer storage
======= =============== ===

Observation input
~~~~~~~~~~~~~~~~~

Format:

::

  {.dbc filename} {time resolution}

The observation data are pulled from a PostgreSQL database. The database connection filename can include a path. The file should provide three queries in the following order:

 1. A query that returns the link_id where observation (gages) are available with the following schema ```CREATE TABLE (link_id integer)```.
 2. A query that returns observation for a given time frame (where begin and end time stamp are parameter) with the following schema ```CREATE TABLE (link_id integer, datetime as integer, discharge real)```.
 3. A query that returns the distance to the border of the domain for the gages with the following schema ```CREATE TABLE (link_id integer, distance real)```.

The time resolution is a floating point number with units in minutes.

Assimilation Window
~~~~~~~~~~~~~~~~~~~

Format:

::

  {num observations}

The duration of the assimilation window expressed in number of time steps.


Forecaster
----------

Running a forescaster with data assimilation requires to run a background simulation with ```asynch``` followed by the analysis with ```assim```. And then generate the forecast using the analysed state as intial conditions.

ON forcing change

RUN asynch for

Notes
-----

.. note::

  The author of these docs is not the primary author of the code so some things may have been lost in translation.

Data assimilation is implemented only for the :ref:`Top Layer Hydrological Model` (``254``). Implementing Data Assimilation requires the user to provide additional model's equations. A more generic method could be used (Jacobian approximation) but would probably be less efficient.

Data assimilation only works with discharge observations (or whatever the first state variable is). This is currently hardcoded but could be extended to support other types of observation such as soil moisture.

Observations should be interpolated to get a better assimilated states (especially for locations that are close to observations).

The larger the assimilation window, the larger is the domain of influence upstream the gages and the better the corrected state. A short assimilation window would only make correction to the links close to the gage and that could induce some ossilations. In Iowa 12 hours, seems to be the sweet spot between computation time and correction.

The solution of the equations at a link depends on the upstreams links and not only the direct parent links. This difference between the forward model and the assimilation model makes Asynch less suitable for solving the system of equations. To be more specific, the partionning of the domain between processors is more senstive since a bad partionning may results in a lot a transfers between procs. Eventually a solver like `CVODES (Sundials) <https://computation.llnl.gov/projects/sundials/cvodes>`_ may be more appropriate.

For small watersheds (N <= 15K links, i.e. Turkey River), ```assim``` works best using serial execution (num procs = 1).

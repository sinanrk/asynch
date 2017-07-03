Data Assimilation
=================

In this section, we discuss the current implementation of the Data Assimilation and the usage of the `assim` CLI.

Introduction
------------

The idea is to combine our hydrological model with stream flow observations. The model used for these notes is the IFC "Golden Standard" model
(252). Section 2 explains the general problem setup. The model and related notation is summarized in Section 3. The mathematical description of combining model outputs with streamflow observations is given in Section 4.1.
Section 4.2 discusses the computations needed for system state sensitivities.
The description is general enough to include a variety of models, though a specic model is used to clarify the presentation.


Building
--------

Data Assimilation is available only if asynch is built with the ``PETSc`` library. Refer to the :ref:`Installation` for more information.




Usage
-----




```
assim turkey_river.gbl turkey_river.das
```



Limitation
----------



### ``NANUK025-ILBOXE53.sh`` is the master script that handles the simulation:

-  creates the different run directories, production, save, restarts, etc
-  fetches and copies all input files and executables into the production directories
-  archives restarts, output files, logs, etc, in the proper directories
-  configure (parses) the template namelists and config files provided here for ``nextsim``,  ``oasis``,  ``opa``, and   ``xios-2.0``, and copy them in the production directory
-  finally generates the batch submission script for the simulation segment to go, example of such a script to be found here in the current directory: 

``run_NANUK025-ILBOXE53_20020101_20020131_00245473-00248448__opa108_nxs112.sub``

As soon as ``NANUK025-ILBOXE53.sh`` works well, it is possible to launch it in the background as follows:

​             `$ nohup ./NANUK025-ILBOXE53.sh &` 

it will run in the background and will wait for the current simulation segment to have completed successfully and will then keep on resubmitting new simulation segments... 

In order to work and find all the **template** config files and namelists,``NANUK025-ILBOXE53.sh`` expects the following subdirectory/file structure into the `Namelists` directory to look like this:

    ./Namelists/iodef.xml 
    ./Namelists/namcouple

    ./Namelists/nxs/cpl_run.cfg

    ./Namelists/opa/namelist_cfg 
    ./Namelists/opa/namelist_ref 
    ./Namelists/opa/domain_def_opa.xml 
    ./Namelists/opa/field_def_opa.xml 
    ./Namelists/opa/file_def_opa.xml 

``lb_functions.sh`` is sourced by ``NANUK025-ILBOXE53.sh`` and contains some useful functions (bash)...



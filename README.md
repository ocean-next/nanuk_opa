# nanuk_opa
Sources, namelists, xml files for the OPA setup of NANUK, based on version 3.6 of NEMOGCM


## Structure

```sources/```: contains only NEMO F90 source files

-  ```sources/NEMOGCM_3.6r7088``` contains the reference NEMOGCM distribution. There is no reason to modify anything in here! This is supposed to be frozen forever!!!
So why is it there then? Well, just for us to make sure we are really using the same background NEMO distro. It's easy to make a "diff" between this directory and the one you are using! 

- ```sources/MY_SRC``` this is where all NEMO sources that differ from what is found into ```sources/NEMOGCM_3.6r7088``` should reside. So there are here all the specific changes specific to Claude/Camille old CREG4 setup + all NANUK-specific tweaks...


```control/```: contains namelists and xml files for OPA and XIOS and OASIS.

-  ```control/opa``` contains namelists for OPA

- ```control/xios-2.0``` contains xml XIOS files for ocean fields output control

- ```control/oasis``` contains the "namcouple" OASIS namelist


```compile/```: contains relevant files needed to compile OPA.

- ```compile/ARCH``` contains relevant "architecture files" for compilation on different HPC systems.

- ```compile/cpp_\<CONFIG\>_OPA.fcm``` contains the CPP keys for a given (domain) configuration
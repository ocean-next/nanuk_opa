#!/bin/bash
#############################################################################
#PBS -N ILBOXE53
#PBS -q mpi_8
#PBS -l select=8:ncpus=28:mpiprocs=28:mem=65G
#PBS -l walltime=05:55:00
#PBS -o out_NANUK025-ILBOXE53_20020101_20020131_00245473-00248448__opa108_nxs112.out
#PBS -e err_NANUK025-ILBOXE53_20020101_20020131_00245473-00248448__opa108_nxs112.err
#PBS -m n
#############################################################################

ulimit -s unlimited

# Name of the config file for neXtSIM
CONFIG=cpl_run.cfg

# Reserved memory for the solver (neXtSIM)
MUMPS_MEM=2048

source /usr/share/Modules/3.2.10/init/sh
export MODULEPATH=/usr/share/modules/modulefiles:/appli/modulefiles:/home3/datawork/eolason/modules:/home3/datawork/eolason/modules:/home3/datawork/eolason/modules
module purge
module load neXtSIM.gcc

# For me (I copied in my HOME so need to override neXtSIM directories:
export NEXTSIMDIR=/home3/datahome/lbrodeau/nextsim
export NEXTSIM_DATA_DIR=/home3/datahome/lbrodeau/nextsim/data/data_links
export NEXTSIM_MESH_DIR=/home3/datahome/lbrodeau/nextsim/mesh

export LD_LIBRARY_PATH=${NEXTSIMDIR}/lib:/home3/datawork/eolason/packages/gmsh/lib:/appli/mpt/2.17-p11472/lib:/appli/netCDF/netcdf-4.4.1.1__gcc-6.3.0__nop/lib:/appli/gcc/6.3.0/6.3.0/lib:/appli/gcc/6.3.0/6.3.0/lib64:/home3/datahome/lbrodeau/local/lib:/appli/netCDF/netcdf-4.4.1.1__gcc-6.3.0__nop//lib

echo
echo "`env | grep -i PBS`" > /home3/datahome/lbrodeau/DEV/nemo_conf_manager/TEST_RUN/CREG025/TEST_CREG025_NEXTSIM-OA3_3.6/ILBOXE53/env_PBS.out

cat ${PBS_NODEFILE} > /home3/datahome/lbrodeau/DEV/nemo_conf_manager/TEST_RUN/CREG025/TEST_CREG025_NEXTSIM-OA3_3.6/ILBOXE53/nodefile_job_${PBS_JOBID}.tmp

nb_cores_pbs=`cat ${PBS_NODEFILE} | wc -l`
echo " *** PBS_NODEFILE contains a list of ${nb_cores_pbs} cores!"

list_nodes=`cat ${PBS_NODEFILE} | awk -v n=28 'NR%n==1'`
#
#list_nodes_c=`scontrol show hostname  | paste -d, -s`
#echo ${list_nodes_c} > /home3/datahome/lbrodeau/DEV/nemo_conf_manager/TEST_RUN/CREG025/TEST_CREG025_NEXTSIM-OA3_3.6/ILBOXE53/node_list.out
#echo
#list_nodes=`echo ${list_nodes_c} | sed -e s/','/' '/g`
#
echo
echo "  *** JOB ID => ${PBS_JOBID} "
echo "  *** Nodes to be booked: ${list_nodes} !"

nb_nodes=`echo ${list_nodes} | wc -w`
if [ ! ${nb_nodes} -eq 8 ]; then echo "ERROR: nb_nodes /= NNODES_TOT !!!"; exit; fi

cd /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/

echo; echo "In `pwd`"; echo; \ls -l; echo

echo "######################################################"
echo " *** Nb cores for xios: 4"
echo " ***      Nodes to use: ${list_nodes_xio}"
echo
echo " *** Nb cores for OPA: 108"
echo " ***      JPNI: 10"
echo " ***      JPNJ: 14"
echo " ***      Nodes to use: ${list_nodes_opa}"
echo
echo " *** Nb cores for NXS: 112"
echo " ***      Nodes to use: ${list_nodes_sas}"
echo
echo " ***       DT_OPA: 900"
echo " *** JOB walltime: 05:55:00"
echo "######################################################"
echo
#
rm -f ocean.output ; # so it does not screw the AAAAAA test later on...
#
${MPI_LAUNCH} \
    -n 27 ./opa.exe : \
    -n 1 ./xios_server.exe : \
    -n 27 ./opa.exe : \
    -n 1 ./xios_server.exe : \
    -n 27 ./opa.exe : \
    -n 1 ./xios_server.exe : \
    -n 27 ./opa.exe : \
    -n 1 ./xios_server.exe : \
    -n 28 ./nextsim.exec --config-files=${CONFIG} -mat_mumps_icntl_23 ${MUMPS_MEM} 1>nxs_log_1.out 2>nxs_log_1.err : \
    -n 28 ./nextsim.exec --config-files=${CONFIG} -mat_mumps_icntl_23 ${MUMPS_MEM} 1>nxs_log_2.out 2>nxs_log_2.err : \
    -n 28 ./nextsim.exec --config-files=${CONFIG} -mat_mumps_icntl_23 ${MUMPS_MEM} 1>nxs_log_3.out 2>nxs_log_3.err : \
    -n 28 ./nextsim.exec --config-files=${CONFIG} -mat_mumps_icntl_23 ${MUMPS_MEM} 1>nxs_log_4.out 2>nxs_log_4.err
#
na=`cat /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/ocean.output 2>/dev/null | grep AAAAAAAA | wc -l`

if [ ${na} -eq 3 ]; then
    # => the job has terminated properly!
    echo "20020131"         > /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/0_last_success_date.info
    echo "85"          > /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/0_last_success_jsub.info
    echo "2976"        > /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/0_last_success_nstock.info
    cat /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/time.step > /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/0_last_success_icpt.info
    # Launching rebuilding of OPA files...
    if true; then
        cd /home3/datawork/lbrodeau/NANUK025/NANUK025-ILBOXE53-S/opa/00245473-00248448/
        rebuild_creg025_prod.sh NANUK025-ILBOXE53
    fi
    #
    # Saving OASIS restarts:
    dirr=/home3/datawork/lbrodeau/NANUK025/NANUK025-ILBOXE53-R/oasis/00248448 ; mkdir -p ${dirr}
    mkdir -p /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/bak
    roa_o=ocean.nc ; roa_i=ice.nc
    for rs in "${roa_o}" "${roa_i}"; do
        rsync -avP /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/${rs} ${dirr}/
        rsync -avP /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/${rs} ${dirr}/
        mv -f /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/${rs} /home3/scratch/lbrodeau/tmp/NANUK025/NANUK025-ILBOXE53_prod/bak/00248448_${rs}
    done
fi
sleep 1
exit

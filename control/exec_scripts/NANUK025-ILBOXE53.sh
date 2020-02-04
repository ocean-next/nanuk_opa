#!/bin/bash

source /usr/share/Modules/3.2.10/init/sh
export MODULEPATH=${MODULEPATH}:/home3/datawork/eolason/modules
module purge
module load neXtSIM.gcc

# For me (I copied in my HOME so need to override neXtSIM directories:
#export NEXTSIMDIR=/home3/datawork/eolason/src/nextsim
export NEXTSIMDIR=/home3/datahome/lbrodeau/nextsim
#export NEXTSIMDIR=/home3/datahome/lbrodeau/nextsim_era5
#
export NEXTSIM_DATA_DIR=${NEXTSIMDIR}/data/data_links
export NEXTSIM_MESH_DIR=${NEXTSIMDIR}/mesh

FREQ_RESTARTS=1 ; # frequency of restarts in months

i_skip_link_restarts=0

. ./lb_functions.sh

ltide=false

lmenage_node=false

# Arch with which NEMO was compiled:
ARCH="DATARMOR_gcc_OA3"

IDEBUG=0

NEMOv="3.6r7088"

CID=""

BDY="ORCA025.L75-GJM189-CREG025.L75-BDY"
BDYI=""

CONF="CREG025" ; CONFO="NANUK025"


DT_NXS=450 ; # time-step for neXtSIM, OPA time step must be a multiple of it, coupling time step a multiple of OPA time step !

list_freq_out="1h 1d 5d" ; # time-tag of files to save/archive

# Datarmor:
NBCPN=28 ; # number of cores per node =  MACHINE SPECIFIC !!!
DIR_STOR_READ_ROOT=/home3/datawork/lbrodeau
DIR_STOR_WRIT_ROOT=/home3/scratch/lbrodeau
DIR_STOR_SAVE_ROOT=/home3/datawork/lbrodeau
#
POST_TRTMT="rebuild_creg025_prod.sh"
launch_post_trmt=true
#
#CP="rsync -LavP" ; # a secure way to copy must follow links!
CP="rsync -LvPI" ; # a secure way to copy must follow links! => gives present timestamp, so stays longer on /scratch !!!

# B A T H Y   T O   U S E :
FBATHY="${DIR_STOR_READ_ROOT}/CREG025/CREG025-I/CREG025_bathy_etopo1_gebco1_smoothed_coast_corrected_bering_may11_time_ct20181224.nc"

# ROOT Directory installed by nemo_conf_manager and in which NEMO should be compiled:
DIR_DEP=${HOME}/NEMO/NEMOv${NEMOv}_OA3-neXtSIM${CID}

# Trying to find the ARCH file:
ARCH_FILE=${DIR_DEP}/ARCH/arch-${ARCH}.fcm
echo ; echo " Archi file: ${ARCH_FILE}"; echo
if [ ! -f ${ARCH_FILE} ]; then
    echo "ERROR: you must specify ARCH_FILE => valid path to the ARCH file used to compile NEMO!"
    exit
fi

MPI_HOME=`cat  ${ARCH_FILE} | grep '^%MPI_HOME'  | sed -e "s|%MPI_HOME||g"  -e "s| ||g" | sed -e 's|${HOME}|/home/users/brodeau|g'`
NCDF_HOME=`cat ${ARCH_FILE} | grep '^%NCDF_HOME' | sed -e "s|%NCDF_HOME||g" -e "s| ||g" | sed -e 's|${HOME}|/home/users/brodeau|g'`
HDF5_HOME=`cat ${ARCH_FILE} | grep '^%HDF5_HOME' | sed -e "s|%HDF5_HOME||g" -e "s| ||g" | sed -e 's|${HOME}|/home/users/brodeau|g'`
XIOS_HOME=`cat ${ARCH_FILE} | grep '^%XIOS_HOME' | sed -e "s|%XIOS_HOME||g" -e "s| ||g" | sed -e 's|${HOME}|/home/users/brodeau|g'`
echo "XIOS_HOME=${XIOS_HOME}" ; echo ; sleep 1
OASIS_HOME=`cat ${ARCH_FILE}| grep '^%OASIS_HOME'| sed -e "s|%OASIS_HOME||g" -e "s| ||g"| sed -e 's|${HOME}|/home/users/brodeau|g'`

OPA_EXE="${DIR_DEP}/CONFIG/${CONF}_OPA/BLD/bin/opa.exe"            ; # OPA executable
XIO_EXE="${XIOS_HOME}/bin/xios_server.exe" ; # XIOS executable
NXS_EXE="${NEXTSIMDIR}/model/bin/nextsim.exec"  ; # neXtSIM executable
echo "Executables to use ="
echo " -> ${OPA_EXE}"
echo " -> ${XIO_EXE}"
echo " -> ${NXS_EXE}"
echo


i_copy_forcing_to_scratch=0

if [ ${IDEBUG} -eq 1 ]; then
    #
    ## We use specific fully dedicated node(s) for xios procs if NNODES_XIO > 0:
    #
    # 2 OPA nodes with 1 xios on each:
    #NCORES_OPA=54; NXO=6 ; NYO=10; NPROC_OPA_P_NODE=27
    #NNODES_XIO=0 ;                 NPROC_XIO_P_NODE=1
    #NCORES_NXS=28 ; QUEUE="mpi_3"; NPROC_NXS_P_NODE=28
    #
    # 3 OPA nodes with 1 xios on each:
    #NCORES_OPA=81; NXO=10 ; NYO=10; NPROC_OPA_P_NODE=27
    #NNODES_XIO=0 ;                  NPROC_XIO_P_NODE=1
    #NCORES_NXS=28 ; QUEUE="mpi_4";  NPROC_NXS_P_NODE=28
    #
    # 4 OPA nodes with 1 xios on each:
    NCORES_OPA=108; NXO=10 ; NYO=14; NPROC_OPA_P_NODE=27
    NNODES_XIO=0 ;                   NPROC_XIO_P_NODE=1
    NCORES_NXS=84 ; QUEUE="mpi_7";   NPROC_NXS_P_NODE=28
    #
    # 5 OPA nodes with 1 xios on each:
    #NCORES_OPA=135; NXO=12 ; NYO=15; NPROC_OPA_P_NODE=27
    #NNODES_XIO=0 ;                   NPROC_XIO_P_NODE=1
    ##NCORES_NXS=28 ; QUEUE="mpi_6";   NPROC_NXS_P_NODE=28
    #NCORES_NXS=56 ; QUEUE="mpi_7";   NPROC_NXS_P_NODE=28
    #
    # 6 OPA nodes with 1 xios on each:
    #NCORES_OPA=162; NXO=18 ; NYO=12; NPROC_OPA_P_NODE=27
    #NNODES_XIO=0 ;                   NPROC_XIO_P_NODE=1
    #NCORES_NXS=28 ; QUEUE="mpi_7";   NPROC_NXS_P_NODE=28
    #
    #
    #
    #
    TJOB="00:55:00" ; # Max wall length (in minutes) for one year of simulation
    #
else
    # PRODUCTION:
    # 3 OPA nodes with 1 xios on each:
    #NCORES_OPA=81; NXO=10 ; NYO=10; NPROC_OPA_P_NODE=27
    #NNODES_XIO=0 ;                  NPROC_XIO_P_NODE=1
    #NCORES_NXS=28 ; QUEUE="mpi_4";  NPROC_NXS_P_NODE=28
    #
    # 4 OPA nodes with 1 xios on each:
    NCORES_OPA=108; NXO=10 ; NYO=14; NPROC_OPA_P_NODE=27
    NNODES_XIO=0 ;                   NPROC_XIO_P_NODE=1
    NCORES_NXS=112 ; QUEUE="mpi_8";  NPROC_NXS_P_NODE=28
    #
    TJOB="05:55:00" # Max wall length (in minutes) for one year of simulation
    #
fi #if [ ${IDEBUG} -eq 1 ]

ca=`printf "%03d" ${NCORES_OPA}` ; cb=`printf "%03d" ${NCORES_NXS}`
cstr_info="_opa${ca}_nxs${cb}"

# Forcing first and last years...
Y1=1995 ; Y2=2020
# ---------------
JD=1  ; # current (cummulated) day since begining of simulation
JDy=1 ; # current day of year

SCRIPT_NAME=`basename $0`

echo

# Getting CONFCASE from the name of the script...
CONFCASE=`basename ${0} | sed -e s/'.sh'/''/g`
#CONFR=`echo ${CONFCASE} | cut --delimiter="-" -f 1,1`
#CONF=`echo ${CONFR} | cut --delimiter="_" -f 1,1`
CASE=`echo ${CONFCASE} | cut --delimiter="-" -f 2,2`
#CONFCASE="${CONFO}-${CASE}"
echo " CONFCASE = ${CONFCASE}"
#echo " CONFR = ${CONFR}"
echo " CONF = ${CONF}"
echo " CONFO = ${CONFO}"
echo " CASE = ${CASE}"
echo " CONFCASE = ${CONFCASE}"
echo

HERE=`pwd`


DATA_CONF_DIR=${DIR_STOR_READ_ROOT}/${CONF}/${CONF}-I
echo " DATA_CONF_DIR => ${DATA_CONF_DIR} "


##############################################################################
# Initializing a run from a given restart files:
#  ( not a natural restart from previous submission )
# ----------------------------------------------------------
init_from_rstrt=0 ; # set to something else than 0 if you want to use external
#                   # restart files, at time step "init_from_rstrt" and specify F_R_INI:
l_respect_rstrt_time=true # if set to false and "init_from_rstrt =/ 0" then will treat the restart as
#                          # an initial state to start from and will ignore the date / time consistency
#                          # and will start the 1st of January !!!
cts_rsrt=`printf "%08d" ${init_from_rstrt}`
F_R_INI=/gpfs/scratch/pr1egh00/pr1egh01/eNATL4/eNATL4-BLB401-R/00000360/eNATL4-BLB401_00000360
##############################################################################






echo; echo "About to launch run ${CONFCASE} forced."; echo; echo





TMPDIR=${DIR_STOR_WRIT_ROOT}/tmp/${CONFO}/${CONFCASE}_prod      ; mkdir -p ${TMPDIR}
RSTDIR=${DIR_STOR_SAVE_ROOT}/${CONFO}/${CONFCASE}-R/opa         ; mkdir -p ${RSTDIR}
RSTDIR_NXS=${DIR_STOR_SAVE_ROOT}/${CONFO}/${CONFCASE}-R/nextsim ; mkdir -p ${RSTDIR_NXS}/out
RSTDIR_OA3=${DIR_STOR_SAVE_ROOT}/${CONFO}/${CONFCASE}-R/oasis   ; mkdir -p ${RSTDIR_OA3}
SAVDIR=${DIR_STOR_SAVE_ROOT}/${CONFO}/${CONFCASE}-S/opa         ; mkdir -p ${SAVDIR}
SAVDIR_NXS=${DIR_STOR_SAVE_ROOT}/${CONFO}/${CONFCASE}-S/nextsim ; mkdir -p ${SAVDIR_NXS}
LOGDIR=${SAVDIR}/LOGS ; mkdir -p ${LOGDIR}


cd ${SAVDIR_NXS}
ln -sf ${RSTDIR_NXS}/out restart

cd ${HERE}/

if [ $(( NCORES_OPA + NCORES_NXS)) -le ${NBCPN} ]; then
    NNODES_OPA=1
    NNODES_NXS=0 # we had to chose on of the 2 to be 0!!! we want the sum to be 1 node!
else

    # Number of nodes for OPA:
    NNODES_OPA=$(( NCORES_OPA/NPROC_OPA_P_NODE ))
    if [ ! $(( NCORES_OPA%NPROC_OPA_P_NODE )) -eq 0 ]; then echo; echo "ERROR: NCORES_OPA is not a multiple of NPROC_OPA_P_NODE !!!"; exit; fi

    # Number of nodes for NXS:
    if [ ! $(( NCORES_NXS%NPROC_NXS_P_NODE )) -eq 0 ]; then echo; echo "ERROR: NCORES_NXS is not a multiple of NPROC_NXS_P_NODE !!!"; exit; fi
    NNODES_NXS=$(( NCORES_NXS/NPROC_NXS_P_NODE ))

fi


# How many XIOS processes:
if [ ${NNODES_XIO} -eq 0 ]; then
    # No specific nodes for XIOS, we put xios processes on the nodes for OPA and NXS!
    #   - 1 check if "NPROC_XIO_P_NODE" compatible with what has been booked for OPA and NXS.
    nb_left=$(( NNODES_OPA*NBCPN - NCORES_OPA ))
    NCORES_XIO=$(( NPROC_XIO_P_NODE*NNODES_OPA ))
    if [ ${NCORES_XIO} -gt ${nb_left} ]; then
        echo "ERROR: NPROC_XIO_P_NODE is not compatible with number of booked nodes for NXS+NEMO !!!"
        echo "     => nb. cores expected for XIOS: ${NCORES_XIO} | left from NEMO+NXS: ${nb_left}"
        exit
    fi
else
    NCORES_XIO=$(($((${NNODES_XIO}>1?${NNODES_XIO}:1))*${NPROC_XIO_P_NODE})) ; # $((${NNODES_XIO}*${NPROC_XIO_P_NODE})) or NPROC_XIO_P_NODE if NNODES_XIO==1...
fi

# Total number of cores and nodes required:
NCORES_TOT=$(( NCORES_OPA + NCORES_NXS + NCORES_XIO ))
NNODES_TOT=$(( NNODES_OPA + NNODES_NXS + NNODES_XIO ))

echo " We are going to book ${NNODES_TOT} nodes and use ${NCORES_TOT} cores! ($(( NNODES_TOT*NBCPN - NCORES_TOT )) procs idle)"
cnxo=`printf "%02d" ${NXO}` ; cnyo=`printf "%02d" ${NYO}`
cnxs=`printf "%02d" ${NXS}` ; cnys=`printf "%02d" ${NYS}`
echo "    OPA => ${NNODES_OPA} nodes (${NCORES_OPA} cores => ${cnxo}x${cnyo})"
echo "        => ${NPROC_OPA_P_NODE} OPA processes on each node fully dedicated to OPA!"
echo "    NXS => ${NNODES_NXS} nodes (${NCORES_NXS} cores => ${cnxs}x${cnys})"
echo "        => ${NPROC_NXS_P_NODE} NXS processes on each node fully dedicated to NXS!"
echo "   XIOS => ${NNODES_XIO} nodes (${NCORES_XIO} cores)"
if [ ${NNODES_XIO} -gt 0 ]; then
    echo "        => ${NPROC_XIO_P_NODE} XIOS processes on each node fully dedicated to XIOS!"
else
    echo "        => ${NCORES_XIO} XIOS processes spread on the same node as OPA and NXS..."
fi
echo


RSTRT="false"; IRCTL=2; TSD_INIT="true"; # default, do not touch

# Namelist in the current directory
fnamelist_o_cfg=Namelists/opa/namelist_cfg
fnamelist_o_ref=Namelists/opa/namelist_ref

fnamelist_x_cfg=Namelists/nxs/cpl_run.cfg


# Grabbing some time stepping/restarting information from the namelist
# ********************************************************************

# Getting NSTOCK from the namelist:
#if [ ! -f ${fnamelist_o_cfg} ]; then echo "Where the hell is the namelist??? (${fnamelist_o_cfg})"; exit; fi
#va=( `cat ${fnamelist_o_cfg} | grep nn_stock | grep = | grep -v stockfl | grep -v '^!'` ); NSTOCK=${va[2]}
#echo; echo "NSTOCK = ${NSTOCK}"; echo

# Getting OPA time step from the namelist:
if [ ! -f ${fnamelist_o_cfg} ]; then echo "Where the hell is the namelist??? (${fnamelist_o_cfg})"; exit; fi
va=( `cat ${fnamelist_o_cfg} | grep 'rn_rdt ' | grep = | grep -v '^!'` ); RDT=${va[2]}
DT_OPA=`echo ${RDT} | sed -e s/'\.'//g`; echo

if [ $(( DT_OPA % DT_NXS )) -ne 0 ]; then
    echo " PROBLEM: OPA time step (${DT_OPA}) is not a multiple of neXtSIM time step (${DT_NXS}) !!!"; exit
fi

DT_CPL=${DT_OPA}

echo
echo " ***  OPA     time step = ${DT_OPA}"
echo " ***  neXtSIM time step = ${DT_NXS}"
echo " *** coupling time step = ${DT_CPL}"
echo




JY=`expr ${Y1} + 0` ; nbdy=`nb_day_in_year ${JY}`

# How many time steps in 1 day and 1 year:
NDT1D=$(((24*3600)/${DT_OPA})) ; echo " => ${NDT1D} time steps per day"
NDT1Y=$((${nbdy}*24*3600/${DT_OPA})) ; echo " => ${NDT1Y} time steps per year"


echo; echo " *** Restart frequency fixed to ${FREQ_RESTARTS} months!"; echo
if   [ ${FREQ_RESTARTS} -eq 6 ]; then
    VDTRSTR=( "0531" )
elif [ ${FREQ_RESTARTS} -eq 1 ]; then
    VDTRSTR=( "0131" "0228" "0331" "0430" "0531" "0630" "0731" "0831" "0930" "1031" "1130" "1231" )
    if [ `lb_is_leap ${JY}` -eq 1 ]; then
        VDTRSTR=( "0131" "0229" "0331" "0430" "0531" "0630" "0731" "0831" "0930" "1031" "1130" "1231" )
    fi
fi
nbrpy=`echo ${VDTRSTR[*]} | wc -w`  ; # number of restarts per year    
echo "  ==> Number of restarts per year: nbrpy=${nbrpy} !"



# Getting neXtSIM time step from the namelist:
if [ ! -f ${fnamelist_x_cfg} ]; then echo "Where the hell is `basename ${fnamelist_x_cfg}`???"; exit; fi

if [ ${init_from_rstrt} -eq 0 ]; then
    l_respect_rstrt_time=true ; # we overide any possible value when init_from_rstrt=0
fi

CN_OCERST_INDIR='.' ; # just default so it's never '' ...




#  ###############################
#  L O O P   A L O N G   Y E A R S
#  ###############################

JYM1=${JY}
jsubm1=0
jsub=1
irstrt=0
icpt=0


while [ ${JY} -le ${Y2} ]; do
    
    echo
    lb_wait_dtrmr ${CASE}
    echo
    
    if [ -f ${TMPDIR}/time.step ]; then
        for ffi in "jsub" "icpt" "date" "nstock"; do
            if [ ! -f ${TMPDIR}/0_last_success_${ffi}.info ]; then
                echo " PROBLEM: [time.step] was found but [0_last_success_${ffi}.info] is missing!!!"; exit
            fi
        done
        DATE_b=`cat ${TMPDIR}/0_last_success_date.info | sed -e "s/ //g"`
        echo "Date at the end of last go is: ${DATE_b}"
        JY=`echo ${DATE_b} | cut -c1-4` ; # changing current year accordingly
    fi

    nbdy=`nb_day_in_year ${JY}`
    
    if [ ${FREQ_RESTARTS} -eq 1 ]; then VDTRSTR=( "0131" "0228" "0331" "0430" "0531" "0630" "0731" "0831" "0930" "1031" "1130" "1231" ); fi
    
    echo
    echo "###########################################################################################"
    if [ -f ${TMPDIR}/time.step ]; then
        jsubm1=`cat ${TMPDIR}/0_last_success_jsub.info` ; jsub=$((jsubm1+1))
        icpt=`cat ${TMPDIR}/0_last_success_icpt.info`
        echo " *** From former go (sub #${jsubm1}) we get: jsub = ${jsub} | icpt = ${icpt} !"
    fi

    if [ ${FREQ_RESTARTS} -eq 1 ] && [ `lb_is_leap ${JY}` -eq 1 ]; then VDTRSTR=( "0131" "0229" "0331" "0430" "0531" "0630" "0731" "0831" "0930" "1031" "1130" "1231" ); fi

    # jidx (C-convention) => index to access the next restart date in vactor "VDTRSTR"

    jidx=`expr ${jsubm1} % ${nbrpy}`
    if [ ${jidx} -lt 0 ]; then echo "PROBLEM: jidx < 0 !!! (${jidx})"; exit; fi
    jidx_b=`expr ${jidx} - 1`
    if [ ${jidx} -eq 0 ]; then jidx_b=`expr ${nbrpy} - 1`; fi
    echo " *** jidx , jidx_b = ${jidx} , ${jidx_b}"
    
    inewyear=0

    iday_rstrt_b=`MMDD_to_day_year ${VDTRSTR[${jidx_b}]} ${JY}` ;
    if [ ${jsub} -gt 1 ]; then
        # This is not 1st submission we're gonna use a restart...
        echo " ==> previous restart creation was at day # ${iday_rstrt_b} of year ${JY}! (${JY}${VDTRSTR[${jidx_b}]})"
        JDy=$((iday_rstrt_b+1))
        JD_last=${iday_rstrt_b}; echo "Last day of year completed: ${JD_last}"
        #
        if [ ${JD_last} -eq ${nbdy} ]; then            
            JY=`expr ${JY} + ${jsub} % ${nbrpy}`
            echo " **** Will start new year => ${JY} ! ****"
            JYM1=`expr ${JY} - 1`
            nbdy=`nb_day_in_year ${JY}`
            if [ ${FREQ_RESTARTS} -eq 1 ] && [ `lb_is_leap ${JY}` -eq 1 ]; then VDTRSTR=( "0131" "0229" "0331" "0430" "0531" "0630" "0731" "0831" "0930" "1031" "1130" "1231" ); fi
            inewyear=1
            JDy=1
            iday_rstrt_b=0
        fi
        #
        SDATE="${JY}`day_year_to_MMDD ${JDy} ${JY}`"
        irstrt=1
        init_from_rstrt=0
        l_respect_rstrt_time=true
        #
    elif [ ${jsub} -eq 1 ]; then
        # This is the first submission!
        JDy=1
        SDATE="${JY}0101" # WE ALWAYS START the 1st of JANUARY !!!
        iday_rstrt_b=0
    else
        echo " WTF??? jsub can only >= 1 !"; exit
    fi
    #        
    # SDATEM1 ?
    if [ ${jsub} -eq 2 ]; then SDATEM1="${JY}0101"; fi
    if [ ${jsub} -ge 3 ]; then
        ii=`MMDD_to_day_year ${VDTRSTR[$((jidx_b-1))]} ${JY}`
        SDATEM1="${JY}`day_year_to_MMDD $((ii+1)) ${JY}`"
        if [ ${inewyear} -eq 1 ]; then
            ii=`MMDD_to_day_year ${VDTRSTR[$((jidx_b-1))]} ${JYM1}`
            SDATEM1="${JYM1}`day_year_to_MMDD $((ii+1)) ${JYM1}`"
        fi
    fi
    #
    echo " => this go is going to start at day # ${JDy} of year ${JY} => SDATE = ${SDATE}"
    JDyend=`MMDD_to_day_year ${VDTRSTR[${jidx}]} ${JY}`
    EDATE="${JY}${VDTRSTR[${jidx}]}" ; echo " => this go is going to stop at day # ${JDyend} of year ${JY} => EDATE = ${EDATE}"
    if [ ${jsub} -ge 2 ]; then
        EDATEM1="${JY}${VDTRSTR[${jidx_b}]}"
        if [ ${inewyear} -eq 1 ]; then EDATEM1="${JYM1}${VDTRSTR[${jidx_b}]}"; fi
    fi
    #
    iday_rstrt=`MMDD_to_day_year ${VDTRSTR[${jidx}]} ${JY}` ; echo " ==> next restart is at ${JY}${VDTRSTR[${jidx}]}, this is day # ${iday_rstrt} of year ${JY}!"
    NB_DAYS_TO_GO=$((iday_rstrt-iday_rstrt_b+ifix)) ;  echo " ===> ${NB_DAYS_TO_GO} days of simulation are going to be completed in this go!"
    NSTOCK=$((NB_DAYS_TO_GO*24*3600/DT_OPA))
    if [ $(((NB_DAYS_TO_GO*24*3600)%DT_OPA)) -ne 0 ]; then "ERROR: cannot get a round NSTOCK!!!"; exit; fi
    echo " ====> nn_stock = ${NSTOCK} time steps to go!" ; echo "  SDATEM1 and EDATEM1 = ${SDATEM1}, ${EDATEM1}"
    echo "###########################################################################################"
    echo
    #    exit;#lolo
    cd ${HERE}/

    # First we need to know if this run is a re-submission
    # ====================================================

    if [ ! -f ${TMPDIR}/time.step ]; then
        irstrt=0
        if [ ${init_from_rstrt} -gt 0 ]; then
            CN_OCERST_INDIR=`dirname ${F_R_INI}`
            fbr=`basename ${F_R_INI}`
            if [ -d ${CN_OCERST_INDIR} ]; then
                echo "INITIAL CONDITION: Using restart files found into ${CN_OCERST_INDIR} !"
            else
                echo "PROBLEM: initial restart directory not found ! (${CN_OCERST_INDIR})"; exit
            fi

            CN_OCERST_IN=${fbr}_restart_oce

            list_rest_o=`ls ${CN_OCERST_INDIR}/${CN_OCERST_IN}*.nc`
            if [ "${list_rest_o}" = "" ]; then
                echo "PROBLEM: no restarts found in directory: ${CN_OCERST_INDIR} !"; exit
            fi

            if ${l_respect_rstrt_time}; then
                # Getting date in restart
                ftst=${CN_OCERST_INDIR}/${CN_OCERST_IN}_0000.nc
                if [ ! -f ${ftst} ]; then
                    ftst=${CN_OCERST_INDIR}/${CN_OCERST_IN}.nc
                    if [ ! -f ${ftst} ]; then
                        echo "PROBLEM: found neither ${CN_OCERST_IN}_0000.nc nor ${CN_OCERST_IN}.nc into ${CN_OCERST_INDIR}/ !"; exit
                    fi
                fi
                kt=`ncdump -v kt ${ftst} | grep 'kt = ' | cut -d "=" -f2 | sed -e s/' '/''/g -e s/';'/''/g`
                SDATEm1=`ncdump -v ndastp ${ftst} | grep 'ndastp = ' | cut -d "=" -f2 | sed -e s/' '/''/g -e s/';'/''/g`
                JY=`echo ${SDATEm1} | cut -c1-4` ; JYend=${JY} ; nbdy=`nb_day_in_year ${JY}`
                echo " => restart_oce says last day done was ${SDATEm1} and kt=${kt} !"
                JDYm1=$((${kt}/${NDT1D}))
                echo " => we have virtually completed ${JDYm1} days since ${JY}0101! "
                JDy=$((${JDYm1}+1))
                JDyend=$((${JDy}+${NB_DAYS_TO_GO}-1))
                echo " => so will start at day # ${JDy} of year ${JY}"
                SDATE="${JY}`day_year_to_MMDD $((${JDy}>1?${JDy}:1)) ${JY}`"
                EDATE=`printf "%04d" ${JY}``day_year_to_MMDD $((${JDyend}>1?${JDyend}:1)) ${JY}`
                # test this SDATE is consistent with date read into restart:
                sdt="${JY}`day_year_to_MMDD $((${JDYm1}>1?${JDYm1}:1)) ${JY}`"
                if [ ! "${sdt}" = "${SDATEm1}" ]; then echo "PROBLEM: dates deduced from restart are fucked up!"; exit; fi
                echo " => JDy, SDATE, EDATE deduced from restart => ${JDy}, ${SDATE}, ${EDATE} !"
                echo
            fi


            RSTRT="true"
            TSD_INIT="false"
        fi

    else

        echo ; echo "OK, we found a ${TMPDIR}/time.step !"; echo

        # Abort if the ocean.output looks suspicious:
        na=`cat ${TMPDIR}/ocean.output | grep AAAAAAAA | wc -l`
        if [ ${na} -eq 3 ]; then
            echo " ocean.output looks OK! " ; echo
        else
            echo "PROBLEM: the last ocean.output looks suspicious! (couldn't see the 3 'AAAAAAAA' rows"
            echo "ABORTING!!! (${TMPDIR}/ocean.output) / na = ${na}" ; echo
            exit
        fi

        # Abort if we find the proof of an borted run:
        ee=`\ls ${TMPDIR}/output.abort_* 2>/dev/null`
        if [ ! "${ee}" = "" ]; then
            echo "PROBLEM: found some output.abort_* files in ${TMPDIR}!!!"
            echo "ABORTING!!!"; echo
            exit
        fi

        citend_last_go=`printf "%08d" ${icpt}`
        CN_OCERST_INDIR="${RSTDIR}/${citend_last_go}"
        CN_OCERST_IN="${CONFCASE}_${citend_last_go}_restart_oce"
        #
        rm -f restart_*.nc
        fto=${CN_OCERST_INDIR}/${CN_OCERST_IN}
        i_happy_with_restart=0
        #
        if [ ${i_skip_link_restarts} -eq 1 ]; then
            echo " *** skipping linking of restarts because i_skip_link_restarts=${i_skip_link_restarts} !!!"
            i_happy_with_restart=1
        fi
        #
        while [ ${i_happy_with_restart} -eq 0 ]; do

            if [ -f ${fto}_0000.nc -o -f ${fto}.nc ]; then
                echo "Restart files are into ${CN_OCERST_INDIR}, good!"; echo
                i_happy_with_restart=1
                cd ${TMPDIR}/
                echo
            else
                echo "Restart files are not into ${CN_OCERST_INDIR}/ !"
                if [ -f ${TMPDIR}/`basename ${fto}_0000.nc` ]; then
                    echo " => but they actually are into ${TMPDIR}/ !"
                    echo "    => moving them into ${CN_OCERST_INDIR}/ !"
                    mkdir -p ${CN_OCERST_INDIR}
                    mv ${TMPDIR}/${CN_OCERST_IN}*.nc ${CN_OCERST_INDIR}/
                else
                    echo "Hey! No restart files found in ${TMPDIR}/ either !!!"
                    echo "   => `basename ${fto}_0000.nc` ?"
                    exit
                fi
            fi

        done

        echo; echo
        # neXtSIM restarts, linking from "${RSTDIR_NXS}/out" to "${RSTDIR_NXS}":
        cd ${RSTDIR_NXS}/
        echo "Expected neXtSIM restarts:"
        for cnr in "mesh" "field"; do
            frnxs_b="${RSTDIR_NXS}/out/${cnr}_${SDATE}T000000Z.bin"
            frnxs_d="${RSTDIR_NXS}/out/${cnr}_${SDATE}T000000Z.dat"
            echo "  -> ${frnxs_b}"; echo "  -> ${frnxs_d}"
            if [ ! -f ${frnxs_b} ] || [ ! -f ${frnxs_d} ]; then echo "ERROR: they are missing!!!"; exit; fi
            ln -sf ${frnxs_b} ./${cnr}_final.bin ; ln -sf ${frnxs_d} ./${cnr}_final.dat ; #because "basename=final" into "cpl_run.cfg"...
        done
        echo "  => done with neXtSIM restarts!"; echo; echo
        cd ${TMPDIR}/

    fi ; # if [ ! -f ${TMPDIR}/time.step ]

    echo
    echo " Will start at day ${SDATE} and end at day ${EDATE}"
    echo " => this is submission # ${jsub}!"
    echo


    if [ ${jsub} -gt 1 ]; then
        # Everything is okay for a normal restart, getting some info about
        # time-steps from revious go:
        # => checking backup of namelist from previous go in "SAVDIR", because
        # namelist_cfg in production dir might be already overwritten, so not
        # reliable
        fnm="${SAVDIR}/namelist_cfg.nsub`printf "%04d" $((${jsub}-1))`" ; #!!! this is the parsed namelist used by last run, so no spaces!
        echo
        if [ ! -f ${fnm} ]; then echo "PROBLEM: previous restart was not found! (${fnm})"; exit; fi
        nn=`cat ${fnm} | grep nn_it000 | grep = | grep -v '^!' | cut -d "=" -f2`
        NN_IT000=`printf "%08d" ${nn}`
        nn=`cat ${fnm} | grep nn_itend | grep = | grep -v '^!' | cut -d "=" -f2`
        NN_ITEND=`printf "%08d" ${nn}`
        echo " *** According to previous namelist (${fnm}):"
        echo "    => Time steps completed during previous go:"
        echo "       from NN_IT000=${NN_IT000} to NN_ITEND=${NN_ITEND}"
        cit_prev="${NN_IT000}-${NN_ITEND}"
        echo ${cit_prev} ; echo
    fi


    CJY=`printf "%04d" ${JY}` ; CJYend=`printf "%04d" ${JYend}`



    # Shall we rely on a restart?
    l_start_from_restart=false
    if [ ${irstrt} -eq 1 ] || [ ${init_from_rstrt} -gt 0 ]; then l_start_from_restart=true; fi



    # Time to create the namelist
    # ---------------------------

    if [ ${init_from_rstrt} -gt 0 ] && ${l_respect_rstrt_time} ; then
        icpt=${init_from_rstrt}
    fi
    IT000=$((${icpt}+1)) ; ITEND=$((${icpt}+${NSTOCK}))
    CITEND=`printf "%08d" ${ITEND}`

    if ${l_start_from_restart}; then
        RSTRT="true" ; IRCTL=2 ; TSD_INIT="false"
        if ! ${l_respect_rstrt_time} ; then IRCTL=0 ; fi
        #
        # Importing OASIS restarts:
        echo; echo "Importing last OASIS restarts!"
        dirr=${RSTDIR_OA3}/${citend_last_go}
        echo " => from dir: ${dirr}/"
        if [ ! -d ${dirr} ]; then echo "PROBLEM: no such directory!!!"; exit; fi
        for coar in "ocean.nc" "ice.nc"; do
            if [ ! -f ${dirr}/${coar} ]; then echo "PROBLEM: no ${dirr}/${coar} !!!"; exit; fi
            ${CP}  ${dirr}/${coar} .
            ln -sf ${dirr}/${coar} ${coar}.lnk
        done
        #
    fi

    # OPA will save its restarts into:
    CN_OCERST_OUTDIR="${RSTDIR}/${CITEND}" ; mkdir -p ${CN_OCERST_OUTDIR}

    # OPA
    # -----
    sed -e "s/<JSUB>/${jsub}/g" -e "s/<CONFCASE>/${CONFCASE}/g" -e "s/<IT000>/${IT000}/g" \
        -e "s/<ITEND>/${ITEND}/g" -e "s/<NN_STOCK>/${NSTOCK}/g" \
        -e "s/<NDATE0>/${SDATE}/g" -e "s/<RSTRT>/${RSTRT}/g" -e "s/<IRCTL>/${IRCTL}/g" \
        -e "s|<CN_OCERST_INDIR>|${CN_OCERST_INDIR}|g" -e "s|<CN_OCERST_IN>|${CN_OCERST_IN}|g" \
        -e "s|<CN_OCERST_OUTDIR>|${CN_OCERST_OUTDIR}|g" \
        -e "s/<TSD_INIT>/${TSD_INIT}/g" -e "s/<NEMOv>/${NEMOv}/g" -e "s/<CONF>/${CONF}/g" \
        -e "s/<JPNI>/${NXO}/g" -e "s/<JPNJ>/${NYO}/g" -e "s/<JPNIJ>/${NCORES_OPA}/g" \
        -e "s/<BDY>/${BDY}/g" -e "s/<BDYI>/${BDYI}/g" \
        ${HERE}/${fnamelist_o_cfg} > ${TMPDIR}/namelist_cfg.human

    ${CP} ${HERE}/${fnamelist_o_ref} ${TMPDIR}/namelist_ref

    # Generating the CFG namelist, void of all crap:
    cat ${TMPDIR}/namelist_cfg.human | grep -o '^[^!]*' | sed '/^\s*$/d' | sed -e s/' '/''/g > \
        ${TMPDIR}/namelist_cfg

    fnmo=${TMPDIR}/namelist_cfg
    echo " *** ${fnmo} generated!"; echo

    fnmob="${SAVDIR}/namelist_cfg.nsub`printf "%04d" ${jsub}`"
    echo " => backing it up as: ${fnmob} !"
    ${CP} ${fnmo} ${fnmob}
    echo


    # NXS
    # ---
    TIME_INIT="`echo ${SDATE}|cut -c1-4`-`echo ${SDATE}|cut -c5-6`-`echo ${SDATE}|cut -c7-8`"
    sed -e "s|<SAVDIR_NXS>|${SAVDIR_NXS}|g" -e "s|<RSTDIR_NXS>|${RSTDIR_NXS}|g" \
        -e "s|<DT_NXS>|${DT_NXS}|g" -e "s|<DT_CPL>|${DT_CPL}|g" \
        -e "s|<NB_DAYS_TO_GO>|${NB_DAYS_TO_GO}|g" -e "s|<Y1>|${Y1}|g" \
        -e "s|<TIME_INIT>|${TIME_INIT}|g" -e "s/<RSTRT>/${RSTRT}/g" \
        ${HERE}/${fnamelist_x_cfg} > ${TMPDIR}/`basename ${fnamelist_x_cfg}`

    fnmx=${TMPDIR}/`basename ${fnamelist_x_cfg}`
    echo " *** ${fnmx} generated!"; echo

    fnmxb="${SAVDIR}/`basename ${fnamelist_x_cfg}`.nsub`printf "%04d" ${jsub}`"
    echo " => backing it up as: ${fnmxb} !"
    ${CP} ${fnmx} ${fnmxb}
    echo


    # OASIS
    # -----
    cruntime=`printf "%08d" $((NB_DAYS_TO_GO*3600*24))`
    sed -e "s|<DT_CPL>|${DT_CPL}|g" -e "s|<DT_OPA>|${DT_OPA}|g" -e "s|<DT_NXS>|${DT_NXS}|g" \
        -e "s|<RUNTIME>|${cruntime}|g" \
        ${HERE}/Namelists/namcouple > ${TMPDIR}/namcouple

    

    # Getting the it info from freshly generated namelist:
    nn=`cat ${fnmo} | grep nn_it000 | grep = | grep -v '^!' | cut -d "=" -f2` ; NN_IT000=`printf "%08d" ${nn}`
    nn=`cat ${fnmo} | grep nn_itend | grep = | grep -v '^!' | cut -d "=" -f2` ; NN_ITEND=`printf "%08d" ${nn}`
    echo " *** Time steps that will be completed during go to be launched:"
    echo "       from NN_IT000=${NN_IT000} to NN_ITEND=${NN_ITEND}"
    cit_now="${NN_IT000}-${NN_ITEND}"


    if [ ! "${CITEND}" = "${NN_ITEND}" ]; then echo "PROBLEM: ITEND =/ NN_ITEND !!! (${CITEND} ${NN_ITEND})"; exit; fi

    ITEND_prev=$((${ITEND}-${NSTOCK}))
    if ${l_start_from_restart} && ${l_respect_rstrt_time} && [ ${ITEND_prev} -ne ${icpt} ]; then
        echo "PROBLEM: ITEND_prev and icpt disagree!!! ${ITEND_prev} and ${icpt}"; exit
    fi

    # Time tag information:
    CTI=${SDATE}_${EDATE}_${cit_now}

    # Backuping namelist as a log:
    ${CP} ${TMPDIR}/namelist_cfg ${SAVDIR}/namelist.${CTI}


    # Installing required files if not a restart:
    if [ ${irstrt} -eq 0 ]; then

        if [ ! -f ${OPA_EXE} ]; then echo; echo "Compile OPA first for ${CONFCASE}!"; exit; fi
        if [ ! -f ${NXS_EXE} ]; then echo; echo "Compile NXS first for ${CONFCASE}!"; exit; fi
        if [ ! -f ${XIO_EXE} ]; then echo; echo "Compile XIOS first for ${CONFCASE}!"; echo "${XIO_EXE}"; exit; fi

        # Copying executables to TMPDIR:
        for pexe in ${OPA_EXE} ${NXS_EXE} ${XIO_EXE}; do
            fexe=`basename ${pexe}`
            ${CP} ${pexe}                                  ${TMPDIR}/
            ln -sf ${pexe} ${fexe}.lnk ; mv -f ${fexe}.lnk ${TMPDIR}/
        done

        # OASIS stuff:
        if [ ! -d ${DATA_CONF_DIR}/OASIS ]; then echo "Hey! We need some OASIS restarts!!! (${DATA_CONF_DIR}/OASIS)"; exit ; fi
        for coar in "ocean" "ice"; do
            foar="${DATA_CONF_DIR}/OASIS/opa_nextsim/${coar}.nc"
            if [ ! -f ${foar} ]; then echo "ERROR! ${foar} is missing!"; exit ; fi
            ${CP} ${foar} ${TMPDIR}/
        done

        # XIOS xml files:
        ${CP} ${HERE}/Namelists/*.xml                 ${TMPDIR}/
        ${CP} ${HERE}/Namelists/opa/*.xml             ${TMPDIR}/

        # Installing configuration files:
        cd ${DATA_CONF_DIR}/
        list=`\ls *_${CONF}.nc`
        cd ${TMPDIR}/
        for ff in ${list}; do
            fn=`echo ${ff} | sed -e s/"_${CONF}.nc"/".nc"/g`
            echo " ${DATA_CONF_DIR}/${ff}  => ${fn} (${CP} ...)"
            ${CP} ${DATA_CONF_DIR}/${ff} ${fn}
        done
        #
        if [ ! -f ${DATA_CONF_DIR}/EMPave_old.dat ]; then echo "HEY! We need a EMPave_old.dat!"; exit; fi
        ${CP} ${DATA_CONF_DIR}/EMPave_old.dat .
        rm -f mesh_mask_*.nc*


        echo;
        echo "********************************************************************************"
        echo "*             We use this bathymetry:";         echo "  ${FBATHY}"
        ${CP} ${FBATHY} bathy_meter.nc
        echo "********************************************************************************"
        echo

        # Copying continental runoffs:
        mkdir -p RNF
        cd RNF/
        ln -sf ${DATA_CONF_DIR}/RNF/*.nc .
        cd ../
        echo

    else

        # This is a restarted year
        echo "Skipping installation of files..."; echo
        cd ${TMPDIR}/
        # Removing old shit:
        rm -f *.out *.err *.tmp tmp* >/dev/null

    fi  # ${irstrt} -eq 0


    echo
    echo " Installing BDY files for year ${JY} !"
    mkdir -p ./BDY
    cd ./BDY/
    for cy in  $((JY-1)) ${JY} $((JY+1)) ; do
        echo "ln -sf ${DATA_CONF_DIR}/BDY/*${BDY}*_y${cy}.nc ."
        ln -sf ${DATA_CONF_DIR}/BDY/*${BDY}*_y${cy}.nc .
    done
    cd ../
    echo


    # the file_def.xml must be updated at each resubmission:
    dir_s=${SAVDIR}/${cit_now} ; mkdir -p ${dir_s}
    sed -e "s|{CONFIG}|${CONFO}|g" -e "s|{CASE}|${CASE}|g" -e "s|{NDATE0}|${SDATE}|g" -e "s|{SAVDIR}|${dir_s}|g" \
        ${HERE}/Namelists/opa/file_def_opa.xml > ${TMPDIR}/file_def_opa.xml
    echo "${HERE}/Namelists/opa/file_def_opa.xml > ${TMPDIR}/file_def_opa.xml"; echo

    if [ ${NB_DAYS_TO_GO} -eq 0 ] && [ "${SDATE}" = "${JY}0101" ]; then
        CTIM1=${SDATE}_${EDATE} ; # same as CTI !!!
    else
        # Normal case:
        CTIM1=${SDATEM1}_${EDATEM1}
    fi


    # Time to launch the nemo executable:
    # ===================================

    cd ${HERE}/

    for ff in "out_${CONFCASE}_*.out" "err_${CONFCASE}_*.err" "run_${CONFCASE}_*.tmp" "node_list.out"; do
        mv -f ${ff} ${LOGDIR}/ 2>/dev/null
    done

    # Create the submission script:
    # =============================

    sub_scr="run_${CONFCASE}_${CTI}_${cstr_info}.tmp"
    rm -f ${sub_scr}

    cat > ${sub_scr} <<EOF
#!/bin/bash
#############################################################################
#PBS -N ${CASE}
#PBS -q ${QUEUE}
#PBS -l select=${NNODES_TOT}:ncpus=${NBCPN}:mpiprocs=${NBCPN}:mem=65G
#PBS -l walltime=${TJOB}
#PBS -o out_${CONFCASE}_${CTI}_${cstr_info}.out
#PBS -e err_${CONFCASE}_${CTI}_${cstr_info}.err
#PBS -m n
#############################################################################

ulimit -s unlimited

# Name of the config file for neXtSIM
CONFIG=cpl_run.cfg

# Reserved memory for the solver (neXtSIM)
MUMPS_MEM=2048

source /usr/share/Modules/3.2.10/init/sh
export MODULEPATH=${MODULEPATH}:/home3/datawork/eolason/modules
module purge
module load neXtSIM.gcc

# For me (I copied in my HOME so need to override neXtSIM directories:
export NEXTSIMDIR=${NEXTSIMDIR}
export NEXTSIM_DATA_DIR=${NEXTSIM_DATA_DIR}
export NEXTSIM_MESH_DIR=${NEXTSIM_MESH_DIR}

export LD_LIBRARY_PATH=\${NEXTSIMDIR}/lib:${LD_LIBRARY_PATH}

echo
echo "\`env | grep -i PBS\`" > ${HERE}/env_PBS.out

cat \${PBS_NODEFILE} > ${HERE}/nodefile_job_\${PBS_JOBID}.tmp

nb_cores_pbs=\`cat \${PBS_NODEFILE} | wc -l\`
echo " *** PBS_NODEFILE contains a list of \${nb_cores_pbs} cores!"

list_nodes=\`cat \${PBS_NODEFILE} | awk -v n=${NBCPN} 'NR%n==1'\`
#
#list_nodes_c=\`scontrol show hostname ${SLURM_NODELIST} | paste -d, -s\`
#echo \${list_nodes_c} > ${HERE}/node_list.out
#echo
#list_nodes=\`echo \${list_nodes_c} | sed -e s/','/' '/g\`
#
echo
echo "  *** JOB ID => \${PBS_JOBID} "
echo "  *** Nodes to be booked: \${list_nodes} !"

nb_nodes=\`echo \${list_nodes} | wc -w\`
if [ ! \${nb_nodes} -eq ${NNODES_TOT} ]; then echo "ERROR: nb_nodes /= NNODES_TOT !!!"; exit; fi

cd ${TMPDIR}/

echo; echo "In \`pwd\`"; echo; \ls -l; echo

echo "######################################################"
echo " *** Nb cores for xios: ${NCORES_XIO}"
echo " ***      Nodes to use: \${list_nodes_xio}"
echo
echo " *** Nb cores for OPA: ${NCORES_OPA}"
echo " ***      JPNI: ${NXO}"
echo " ***      JPNJ: ${NYO}"
echo " ***      Nodes to use: \${list_nodes_opa}"
echo
echo " *** Nb cores for NXS: ${NCORES_NXS}"
echo " ***      Nodes to use: \${list_nodes_sas}"
echo
echo " ***       DT_OPA: ${DT_OPA}"
echo " *** JOB walltime: ${TJOB}"
echo "######################################################"
echo
#
rm -f ocean.output ; # so it does not screw the AAAAAA test later on...
#
\${MPI_LAUNCH} \\
EOF
    for ii in $(seq 1 ${NNODES_OPA}); do
        echo "    -n ${NPROC_OPA_P_NODE} ./opa.exe : \\" >> ${sub_scr}
        echo "    -n ${NPROC_XIO_P_NODE} ./xios_server.exe : \\" >> ${sub_scr}
    done
    cec=" : \\"
    for ii in $(seq 1 ${NNODES_NXS}); do
        if [ ${ii} -eq ${NNODES_NXS} ]; then cec=""; fi
        echo "    -n ${NPROC_NXS_P_NODE} ./nextsim.exec --config-files=\${CONFIG} -mat_mumps_icntl_23 \${MUMPS_MEM} 1>nxs_log_${ii}.out 2>nxs_log_${ii}.err${cec}" >> ${sub_scr}
    done
    #
    cat >> ${sub_scr} <<EOF
#
na=\`cat ${TMPDIR}/ocean.output 2>/dev/null | grep AAAAAAAA | wc -l\`

if [ \${na} -eq 3 ]; then
    # => the job has terminated properly!
    echo "${EDATE}"         > ${TMPDIR}/0_last_success_date.info
    echo "${jsub}"          > ${TMPDIR}/0_last_success_jsub.info
    echo "${NSTOCK}"        > ${TMPDIR}/0_last_success_nstock.info
    cat ${TMPDIR}/time.step > ${TMPDIR}/0_last_success_icpt.info
    # Launching rebuilding of OPA files...
    if ${launch_post_trmt}; then
        cd ${dir_s}/
        ${POST_TRTMT} ${CONFCASE}
    fi
    #
    # Saving OASIS restarts:
    dirr=${RSTDIR_OA3}/${CITEND} ; mkdir -p \${dirr}
    mkdir -p ${TMPDIR}/bak
    roa_o=ocean.nc ; roa_i=ice.nc
    for rs in "\${roa_o}" "\${roa_i}"; do
        rsync -avP ${TMPDIR}/\${rs} \${dirr}/
        rsync -avP ${TMPDIR}/\${rs} \${dirr}/
        mv -f ${TMPDIR}/\${rs} ${TMPDIR}/bak/${CITEND}_\${rs}
    done
fi
sleep 1
exit
EOF
    # Script is done!

    echo; echo "Submitting batch job ${sub_scr} ..."
    echo "qsub ${sub_scr}"
    chmod +x ${sub_scr}
    qsub ${sub_scr}

    if [ ${IDEBUG} -eq 1 ]; then echo " We exit because IDEBUG=${IDEBUG}!"; exit; fi

    echo ; echo " Sleeping 15 seconds..."; sleep 15; echo

    # Waiting until the batch job is over...
    lb_wait_dtrmr ${CASE}

    cd ${TMPDIR}/


    echo; echo;

    JYM1=${JY}

done # loop along years...

##############################################################################
#
#      Input file for OASIS3
# 
#     NANUK aka "OPA -- OASIS -- neXtSIM" on CREG025 (source and target) domain
#   
#     2019: E. Olason & L. Brodeau
#
#    Important:
#    ----------
#   <RUNTIME> and <DT_CPL> have to be replaced by the length of the simulation (s) and
#   the coupling frequency (s), respectively.
#
#
###############################################################################
#  
#      Input delimiters have to occupy position 1 to 9 !
#      No blank lines allowed !
#      Length of input lines <= 80 !
#
###############################################################################
#
# NFIELDS : total number of fields being exchanged.
#
 $NFIELDS
         14
 $END
#
 $RUNTIME
           <RUNTIME>
 $END
#
 $NLOGPRT
   0  0
 $END
#
 $STRINGS
#
#
############################################################################
#                      OCEAN  --->>>  ICE
#                      ------------------
###############################################################################
#
# --- start Field 1 --- Sea_surface_temperature [K;K]
#
O_SSTSST I_SST 1 <DT_CPL>  1 ocean.nc EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
AVERAGE
#
# --- end field 1 ---
#
###############################################################################
#
# --- start Field 2 --- sea_surface_salinity  [PSU;PSU]
#
O_SSSal  I_SSS 2 <DT_CPL>  1 ocean.nc EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
AVERAGE
#
# --- end field 2 ---
#
###############################################################################
#
# --- start Field 3 --- zonal_current  [m/s;m/s]
#
O_OCurx1 I_Uocn 2 <DT_CPL>  1 ocean.nc EXPORTED
528 603  528 603 uor2 uor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
AVERAGE
#
#
# --- end field 3 ---
#
###############################################################################
#
# --- start Field 4 --- meridional_current  [m/s;m/s]
#
O_OCury1 I_Vocn 2 <DT_CPL>  1 ocean.nc EXPORTED
528 603  528 603 vor2 vor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
AVERAGE
#
#
# --- end field 4 ---
#
###############################################################################
#
# --- start Field 5 --- SSHeight  [m;m]
#
O_SSHght I_SSH 2 <DT_CPL>  1 ocean.nc EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
AVERAGE
#
#
# --- end field 5 ---
#
###############################################################################
#
# --- start Field 6 --- Fraction of solar net radiation absorbed in the first ocean level
#
O_FraQsr I_FrcQsr 2 <DT_CPL>  1 ocean.nc EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
AVERAGE
#
#
# --- end field 6 ---
#
#
#
#
#
#
#
############################################################################
#                      ICE  --->>>  OCEAN
############################################################################
#
# --- start Field 7 --- eastward  wind stress over ocean at U point
#
I_taux O_OTaux1   23  <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 uor2 uor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
INSTANT
#
# --- end field 7 ---
#
#########################################################################
#
# --- start Field 8 --- northward wind stress over ocean at V point
#
I_tauy O_OTauy1   24  <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 vor2 vor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
INSTANT
#
# --- end field 8 ---
#
############################################################################
#
# --- start Field 9 --- Total E-P -> [ Evaporation - (rain + snow) ]
#
I_fwflux OOEvaMPr 28 <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
ACCUMUL
#
#
# --- end field 9 ---
#
###########################################################################
#
# --- start Field 10 --- Total Non Solar 
#
I_rsnos O_QnsOce 6  <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
ACCUMUL
#
#
# --- end field 10 ---
#
############################################################################
#
# --- start Field 11 --- Total Solar 
#
I_rsso O_QsrOce 7 <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
ACCUMUL
#
# --- end field 11 ---
##########################################################################
#
# --- start Field 12 --- Salt flux
#
I_sfi O_SFLX 7 <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
ACCUMUL
#
# --- end field 12 ---
#
##########################################################################
#
# --- start Field 13 --- Wind stress module       
#
I_taumod  O_TauMod 7 <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
INSTANT
#
# --- end field 13 ---
##########################################################################
#
# --- start Field 14 ---  Sea ice cover      
#
I_sic RIceFrc 7 <DT_CPL>  1  ice.nc  EXPORTED
528 603  528 603 tor2 tor2  SEQ=1 LAG=<DT_CPL>
R  0  R  0
#
LOCTRANS
INSTANT
#
# --- end field 14 ---
##########################################################################
 $END

#!/bin/csh
#
# batobs.sh (updated): read MASS from ../VADAT where the line format is:
#   <value>   [MASS]
# with arbitrary whitespace between the value and [MASS].
#
# Expects: run from MODEL/obs/ so ../VADAT and ../RVTJ exist; $cmfdist and
# $ATOMIC must be set (e.g. in environment or before sourcing/calling ../batch.sh).
#

setenv PROG_CMF_OBS   $cmfdist/exe/cmf_flux.exe

# ***********************************************************************
# Get MASS from ../VADAT
# ***********************************************************************
set vadat = "../VADAT"

if ( ! -e "$vadat" ) then
  echo "ERROR: $vadat not found; cannot read MASS"
  exit 1
endif

# Match any line containing [MASS], print first whitespace-separated field.
# This is robust to arbitrary spaces/tabs between value and [MASS].
set MASS = `awk 'index($0,"[MASS]")>0 {print $1; exit}' "$vadat"`

if ( "$MASS" == "" ) then
  echo "ERROR: Could not parse MASS from $vadat (expected: <value>   [MASS])"
  echo "First few lines of $vadat:"
  head -n 20 "$vadat"
  exit 1
endif

echo "Using MASS = $MASS from $vadat"

# ***********************************************************************
# Perform soft links so Atomic data files can be accessed.
# ***********************************************************************
../batch.sh assign

ln -sf $ATOMIC/misc/strk_list_lemke_3jan01.dat   FULL_STRK_LIST
ln -sf $ATOMIC/misc/lemke_hi.dat                 LEMKE_HI
ln -sf $ATOMIC/misc/hei_ir_strk.dat              HeI_IR_STRK
ln -sf $ATOMIC/misc/ds_strk.dat                  DS_STRK
ln -sf $ATOMIC/misc/bs_hhe.dat                   BS_HHE
ln -sf $ATOMIC/misc/REVISED_LAMBDAS_9apr08       REVISED_LAMBDAS

# ***********************************************************************
# Copy file across with main settings (TRAP_J forced to T).
# ***********************************************************************
sed -e "/\\[TRAP_J\\]/s/.*/    T      [TRAP_J]/" \
  CMF_FLUX_PARAM_INIT > CMF_FLUX_PARAM

rm -rf OUT_FLUX

# ***********************************************************************
# Start the job
# ***********************************************************************
rm -f batobs.log
rm -f error.log

echo "Machine name is :" > 'batobs.log'
uname -a >> 'batobs.log'

echo " " >> 'batobs.log'
echo "Program started on:" >> 'batobs.log'
date >> 'batobs.log'
echo " " >> 'batobs.log'

(exec nice -15 $PROG_CMF_OBS  >>& 'batobs.log')  << END
../RVTJ            [RVTJ]
$MASS              [MASS]
F                  [ONLY_OBS_LINES]
END

echo "Program finished on:" >> 'batobs.log'
date >> 'batobs.log'

# ***********************************************************************
# Clean up
# ***********************************************************************
mv -f OBSFRAME     obs_fin_10
mv -f HYDRO        hydro_fin_10
mv -f MEANOPAC     meanopac_fin
mv -f TIMING       full_timing
mv -f J_COMP       J_COMP_full

rm -f ES_J_CONV
rm -f TRANS_INFO
rm -f fort.*

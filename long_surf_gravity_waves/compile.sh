#! /bin/bash
#set -x

EXE=shallow_water.exe

template=$1

if [[ -z $template ]]; then
    echo 'Error: No template file given.'
    echo usage:
    echo $0 template_file
    exit 1
fi


thisdir=$(pwd)

template=$thisdir/$template

srcdir=$(dirname $thisdir)

FMS_UTILS=$srcdir/fms_shared

FMS_UTILITIES="$FMS_UTILS/include $FMS_UTILS/platform $FMS_UTILS/constants $FMS_UTILS/fms $FMS_UTILS/time_manager $FMS_UTILS/mpp $FMS_UTILS/diag_manager  $FMS_UTILS/memutils FMS_UTILS/constants"

MKMF=$FMS_UTILS/mkmf

mkdir -p $thisdir/exec

cd $thisdir/exec

$MKMF -m Makefile -p $EXE -t $template $thisdir $FMS_UTILITIES  $FMS_UTILS/mpp/include

make "${@:2}"

exit

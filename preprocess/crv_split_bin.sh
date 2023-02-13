# /bin/bash
#
# as long as nothing has changed, then running this script will rerite files from
# the "output" dir to the "temp" dir, then those have to be mv'ed into "upload"
#
# current assumptions are that inputs are in the text format with names
#    RUN_FEB_60_TB202106171229_GeV_Kcnt_BV_900-900-Gain_test_deg_x_z_degF.data
# or of the form:
#    raw.mu2e.CRV_wideband_cosmics.crvaging-003.001053.dat
# and that output files should be like
#    raw.mu2e.CRV_wideband_cosmics.crvaging-001.000042_000.dat
#
#


# where file to be split can be found
export ID=/mnt/d/data
# export ID=/mnt/c/data/process
# where to write splits files
export OT=/mnt/d/data_split
# export OT=/mnt/c/data/splitted
# where to find crv_split.py
export SOURCEDIR=/mnt/d

for FF in $( ls $ID )
#for FF in RUN_FEB_101_TB202201071311_GeV_Kcnt_BV_900-900-900-900-Gain_test_deg_x_z_degF.data
do
    echo "**** $FF"
    INFS=$ID/$FF
    NFIELDS=$( echo $FF | awk -F. '{print NF}' )

    if [ $NFIELDS -eq 6 ]; then  # standard format
        RR=$(echo $FF | awk -F. '{print $5}' | awk -F_ '{print $1}' | sed 's/^0*//' )
        BBB=$(echo $FF | awk -F. '{print $1 "." $2 "." $3 "." $4 }')
        EEE=$(echo $FF | awk -F. '{print $6}')
        FN=$( printf "%s.%06d_000.%s" $BBB $RR $EEE)
    else  # DAQ format
        echo "DAQ format code is not up to date"
        exit 1
        RR=$(echo $FF | awk -F_ '{print $3}' )
        CONF=1
        [[ "$FF" =~ "139" ]] && CONF=2
        FN=$( printf "raw.mu2e.CRV_wideband_cosmics.crvaging-%03d.%06d_000.dat" $CONF $RR )
    fi

    OTFS=$OT/$FN

    # MB
    SZ=$( ls -l $INFS | awk '{ print int( $5/1025/1024 ) }' )
    NE=$( grep -c FEB0 $INFS )
    NL=$( wc -l $INFS | awk '{print $1}' )

    FTYPE="aging"
    NSPLIT=1000
    if [[ "$FF" =~ "crvled" ]]; then
        FTYPE="led"
        NSPLIT=5
    fi
    NOUT=$((NE/NSPLIT))
    [ $NOUT -lt 1 ] && NOUT=1

    echo "Size $SZ lines: $NL Nev: $NE  file_type: $FTYPE nout: $NOUT"

    # split if there are at least two whole potential output files

    if [ $NOUT -le 1 ]; then
        echo "Just copying "
        rm -f $OTFS
        cp $INFS $OTFS
    else
        echo "Splitting ./crv_split_bin.py $FF"

        # note: this code contains file name assumptions

        $SOURCEDIR/crv_split_bin.py $INFS $OTFS $NE $NSPLIT
        #$SOURCEDIR/crv_split_bin.py $INFS $OTFS 1000

    fi


done

exit 0

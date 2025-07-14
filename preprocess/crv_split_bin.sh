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
# export ID=/mnt/d/raw
export ID=/mnt/c/data/raw
# where to write splits files
#export OT=/mnt/d/splitted
export OT=/mnt/c/data/splitted
# where to find crv_split.py
# export SOURCEDIR=/mnt/d
export SOURCEDIR=/mnt/c/data

for FF in $( ls $ID )
#for FF in RUN_FEB_101_TB202201071311_GeV_Kcnt_BV_900-900-900-900-Gain_test_deg_x_z_degF.data
do
    echo "**** $FF"
    INFS=$ID/$FF
    #echo "**** $INFS"
    NFIELDS=$( echo $FF | awk -F. '{print NF}' )
    #echo "**** $NFIELDS"

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
    #echo "**** $OTFS"

    # MB
    SZ=$( ls -l $INFS | awk '{ print int( $5/1025/1024 ) }' )
    #echo "**** $SZ"
    # NE=$( grep -c FEB0 $INFS )
    NE=$( rg -c FEB0 $INFS ) # use ripgrep
    #echo "**** $NE"
    NL=$( wc -l $INFS | awk '{print $1}' )
    #echo "**** $NL"

    FTYPE="aging"
    NSPLITSUPPLIED="${1:-0}"
    NSPLIT=0
    if [ $NSPLITSUPPLIED -le 0 ]; then
        echo "Using default spills per subrun"
        if [[ "$FF" =~ "crvled" ]]; then
            FTYPE="led"
            NSPLITSUPPLIED=5
            echo "Assigning spills per subrun = 5 for file type crvled"
        else
            NSPLITSUPPLIED=250
            echo "Assigning spills per subrun = 250 for file type crvaging"
        fi
    else
        echo "Required spills per subrun: $NSPLITSUPPLIED"
    fi

    NSPLIT=$((NSPLITSUPPLIED))
    NOUT=$((NE/NSPLIT))
    [ $NOUT -lt 1 ] && NOUT=1
    SZPERSUBRUN=$((SZ/NOUT))

    if [ $SZPERSUBRUN -gt 20000 ]; then
        echo "ABORT!!! Subrun files are too huge"
        echo "TOTAL     Size: $SZ    Nspills: $NE"
        echo "PARSED    Size: $SZPERSUBRUN    Nspills: $NSPLIT"
        echo "Please adjust # spills per subrun"
        exit 2
    fi  

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
        #$SOURCEDIR/crv_split_bin.py $INFS $OTFS 100

    fi


done

exit 0

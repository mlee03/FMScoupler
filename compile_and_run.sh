#!/bin/bash

function compile() {

    #1 platform
    #2 target
    #3 compile experiment
    #4 xml
    #5 model
    
    LOG=$5
    
    echo "=========================" > $LOG
    echo "COMPILE `date`" >> $LOG
    echo "COUPLER_HASH $coupler_hash" >> $LOG
    echo "=========================" >> $LOG

    ( fremake -f -x $4 -p $1 -t $2 $3 ) >> $LOG 2>&1
    runscript=$( grep "TO SUBMIT" $LOG | awk '{print $NF}' | tail -n 1 )
    ( bash -c $runscript ) >> $LOG 2>&1
}


function submit() {

    #1 platform
    #2 target
    #3 regression test
    #4 experiment
    #5 xml
    #6 model

    LOG=$6

    echo "=========================" >> $LOG
    echo "SUBMIT $4 $3 `date`" >> $LOG
    echo "COUPLER_HASH $coupler_hash" >> $LOG
    frerun -s --overwrite --no-transfer -x $5 -p $1 -t $2 -r $3 $4 >> $LOG 2>&1
    echo "=========================" >> $LOG

}


# --no-free



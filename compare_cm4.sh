#!/bin/bash

benchmark_dir=/gpfs/f5/gfdl_f/scratch/Mikyung.Lee/CM4/main
test_dir=/gpfs/f5/gfdl_f/scratch/Mikyung.Lee/CM4/$1

expts=( CM4_piControl_C )
prod=ncrc5.intel23-repro-openmp
rest_hist_list=( restart history )

logfile=NCCMP_cm4

echo "===================================" > $logfile
echo `date`                                >> $logfile
echo "testing branch $1"                   >> $logfile
echo "===================================" >> $logfile

for iexpt in ${expts[@]} ; do
    for idir in $(cd $benchmark_dir/$iexpt/$prod/archive && ls -d [0-9]*) ; do
        for rest_hist in restart history ; do
            bdir=$benchmark_dir/$iexpt/$prod/archive/$idir/$rest_hist
            tdir=$test_dir/$iexpt/$prod/archive/$idir/$rest_hist
            echo ""                                           >> $logfile            
            echo "==================================="        >> $logfile
            echo $iexpt $idir $idir $rest_hist                >> $logfile
            echo "==================================="        >> $logfile
            for itar in $bdir/*.tar ; do
                tar -xf $itar -C $bdir
                tar -xf ${itar/$bdir/$tdir} -C $tdir
                for ifile in $bdir/*.nc ; do
                    echo ${ifile/$bdir} >> $logfile
                    nccmp -c 1 -md -f -s  $ifile ${ifile/$bdir/$tdir} >> $logfile 2>&1
                done
                #ls $itar ${itar/$bdir/$tdir} | ardiff -fg     >> $logfile 2>&1 &
            done
        done
    done
done
echo "DONE" >> $logfile

### check ascii
for iexpt in ${expts[@]} ; do
    for idir in $(cd $benchmark_dir/$iexpt/$prod/archive && ls -d [0-9]*) ; do
        bdir=$benchmark_dir/$iexpt/$prod/archive/$idir/ascii
        tdir=$test_dir/$iexpt/$prod/archive/$idir/ascii
        echo ""                                           >> $logfile
        echo "==================================="        >> $logfile
        echo "$iexpt $idir $idir ascii"                   >> $logfile
        echo "==================================="        >> $logfile
        
        for itar in $bdir/*.tar ; do
            tar -xf $itar -C $bdir
            tar -xf ${itar/$bdir/$tdir} -C $tdir
            #check CHECKSUMS in logfile
            ifile=file.$$
            grep "^CHECKSUM" $bdir/*$idir.o* > $ifile
            while read -r iline ; do
                icount=$(grep -c "$iline" $tdir/*$idir.o*)
                if [[ $icount == 0 ]] ; then echo "DIFFER $iline" >> $logfile ; fi
            done < $ifile
            rm $ifile
        done
        
        #check ocean.stats
        for istats in $bdir/*ocean.stats.nc ; do
            echo ${istats/$bdir}                                >> $logfile
            nccmp -c 1 -md -f -s  $istats ${istats/$bdir/$tdir} >> $logfile 2>&1 
        done
        
        #check stocks
        for istocks in $bdir/*stocks.out ; do
            echo ${istocks/$bdir}                >> $logfile
            diff $istocks ${istocks/$bdir/$tdir} >> $logfile 2>&1 
        done
    done
done
echo "DONE" >> $logfile


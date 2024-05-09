#!/bin/bash

source compile_and_run.sh

coupler_hash=$1

platform=ncrc5.intel23-classic
target=repro-openmp
rts=rts

xml_dir=/ncrc/home1/Mikyung.Lee/FMScoupler/tests/xml/AM4/xml
xml=$xml_dir/awg_xanadu.xml

my_xml=${xml}_$coupler_hash
sed "s/REPLACEME/$coupler_hash/g" < $xml > $my_xml

expt_list=( c96L33_am4p0 c96L33_am4p0_cmip6Diag c96L33_am4p0_2010climo
            c96L33_am4p0_1850climo c192L33_am4p0_cmip6Diag )
compile_expt=cm4p12_2023.04

LOG=$PWD/LOG_${coupler_hash}_am4

if [ $# -eq 2 ] ; then
    if [ $2 == 'compile' ] ;  then
        cd $xml_dir
        compile $platform $target $compile_expt $my_xml $LOG
    elif [ $2 == 'run' ] ; then
        cd $xml_dir
        for iexpt in ${expt_list[@]} ; do 
            if [ $? -eq 0 ] ; then submit $platform $target $rts $iexpt $my_xml $LOG  ; fi
        done
    fi
else 
    cd $xml_dir
    compile $platform $target $compile_expt $my_xml $LOG
    if [ $? -eq 0 ] ; then
        cd $xml_dir
        for iexpt in ${expt_list[@]} ; do 
            if [ $? -eq 0 ] ; then submit $platform $target $rts $iexpt $my_xml $LOG ; fi            
        done
    fi
fi    

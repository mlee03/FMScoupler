#!/bin/bash

#set -x 
source compile_and_run.sh

coupler_hash=$1

platform=ncrc5.intel23
target=repro-openmp
rts=rts  #basic

xml_dir=/ncrc/home1/Mikyung.Lee/FMScoupler/tests/xml/CM4_ESM4_OM4/mdt_xml/regression_xmls 
xml=$xml_dir/CM4_piControl_C.xml
new_xml=${xml}_${coupler_hash}

LOG=$PWD/LOG_${coupler_hash}_cm4

sed "s/REPLACEME/$coupler_hash/g" < $xml > $new_xml

compile_expt=cm4_compile
expt=CM4_piControl_C

if [ $# -eq 2 ] ; then
    if [ $2 == 'compile' ] ; then
        cd $xml_dir
        compile $platform $target $compile_expt $new_xml $LOG
    elif [ $2 == 'run' ] ; then
        cd $xml_dir
        submit $platform $target $rts $expt $new_xml $LOG
    fi
else
    cd $xml_dir
    compile $platform $target $compile_expt $new_xml $LOG
    if [ $? -eq 0 ] ; then
        cd $xml_dir
        submit $platform $target $rts $expt $new_xml $LOG
    fi
fi

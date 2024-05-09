#!/bin/bash

#set -x 

source compile_and_run.sh

coupler_hash=$1

platform=ncrc5.intel23-classic
target=repro-openmp
rts=basic_prodsettings

xml_dir=/ncrc/home1/Mikyung.Lee/FMScoupler/tests/xml/SPEAR

xml=$xml_dir/"SPEAR_c192_o1_Control_1850_Q50.xml"
compile_xml=$xml_dir/SPEAR_experiments_compile.xml

my_xml=${xml}_${coupler_hash}
my_compile_xml=${compile_xml}_${coupler_hash}

LOG=$PWD/LOG_${coupler_hash}_spear

sed "s/REPLACEME/$coupler_hash/" < $compile_xml > $my_compile_xml
sed "s/REPLACEME/$coupler_hash/" < $xml > $my_xml

compile_expt=SPEAR_Q2023.04_nonsymMOM6_exec
expt_list=( SPEAR_c192_o1_Control_1850_Q50 )

if [ $# -eq 2 ] ; then
    if [ $2 == 'compile' ] ; then
        cd $xml_dir
        compile $platform $target $compile_expt $my_xml $LOG
    elif [ $2 == 'run' ] ; then
        cd $xml_dir
        for iexpt in ${expt_list[@]} ; do
            if [ $? -eq 0 ] ; then submit $platform $target $rts $iexpt $my_xml $LOG ; fi
        done
    fi
else
    cd $xml_dir
    compile $platform $target $compile_expt $my_xml $LOG
    for iexpt in ${expt_list[@]} ; do
        if [ $? -eq 0 ] ; then submit $platform $target $rts $iexpt $my_xml $LOG ; fi
    done
fi
    

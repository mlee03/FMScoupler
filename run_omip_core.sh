#!/bin/bash

source compile_and_run.sh

coupler_hash=$1

platform=ncrc5.intel23
target=repro
rts=rts

xml_dir=/ncrc/home1/Mikyung.Lee/FMScoupler/tests/xml/CM4_ESM4_OM4/mdt_xml

compile_xml=$xml_dir/awg_include/xml_building_blocks/om4_compile.xml
xml=$xml_dir/regression_xmls/OMIP4_CORE2.xml

my_compile_xml=${compile_xml}_${coupler_hash}
my_xml=${xml}_${coupler_hash}

LOG=$PWD/LOG_${coupler_hash}_omip_core2

sed "s/REPLACEME/$coupler_hash/g" < $compile_xml > $my_compile_xml
sed "s/om4_compile.xml/om4_compile.xml_${coupler_hash}/g" < $xml > $my_xml
sed -i "s/REPLACEME/$coupler_hash/g" $my_xml

expt_list=( OM4p25_IAF_BLING_csf_rerun OM4p5_IAF_BLING_CFC_abio_csf_mle200 )
compile_expt=MOM6_SIS2_compile_FMS2

if [ $# -eq 2 ] ; then
    if [ $2 == 'compile' ] ;  then
        cd $xml_dir/regression_xmls
        compile $platform $target $compile_expt $my_xml $LOG
    elif [ $2 == 'run' ] ; then
        for iexpt in ${expt_list[@]} ; do 
            if [ $? -eq 0 ] ; then
                cd $xml_dir/regression_xmls
                submit $platform $target $rts $iexpt $my_xml $LOG
            fi
        done
    fi
else 
    compile $platform $target $compile_expt $my_xml $LOG
    if [ $? -eq 0 ] ; then
        for iexpt in ${expt_list[@]} ; do 
            cd $xml_dir/regression_xmls
            if [ $? -eq 0 ] ; then
                submit $platform $target $rts $iexpt $my_xml $LOG
            fi
        done
    fi
fi    

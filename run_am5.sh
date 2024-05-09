#!/bin/bash

source compile_and_run.sh

coupler_hash=$1

ver=f5b6r0

platform=ncrc5.intel23-classic
target=repro-openmp
rts=prodrts

xml_dir=/ncrc/home1/Mikyung.Lee/FMScoupler/tests/xml/am5xml

xml=$xml_dir/am5.xml
compile_xml=$xml_dir/xml_experiments/compile.xml

my_xml=${xml}_${coupler_hash}
my_compile_xml=${compile_xml}_${coupler_hash}

LOG=$PWD/LOG_${coupler_hash}_am5

sed "s/REPLACEME/$coupler_hash/g" < $xml > $my_xml
sed "s/REPLACEME/$coupler_hash/g" < $compile_xml > $my_compile_xml

expt_list=( c96L65_am5${ver}_amip c96L65_am5${ver}_pdclim1850F c96L65_am5${ver}_pdclim2010F
            c384L65_am5${ver}_amip c384L65_am5${ver}_pdclim2010F c384L65_am5${ver}_pdclim1850F
            c384L65_am5${ver}_OM4_p25_piControl_noBLING_DynVeg c96L65_am5${ver}_OM4_p25_piControl_noBLING_DynVeg
            c96L65_am5${ver}_amip_cosp )
compile_expt=am5${ver}_compile

if [ $# -eq 2 ] ; then
    if [ $2 == 'compile' ] ;  then
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

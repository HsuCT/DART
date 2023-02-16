#!/bin/bash

let year=2013
let month=03
let day=17
let hour=0
let minute=0
let second=0


while (( hour <= 23)) ; do
str_year=`printf %04d $year`
str_month=`printf %02d $month`
str_day=`printf %02d $day`
str_hour=`printf %02d $hour`
str_minute=`printf %02d $minute`
str_second=`printf %02d $second`


  sed -e "s/YYYY/${str_year}/g"   \
      -e "s/MM/${str_month}/g"    \
      -e "s/DD/${str_day}/g"      \
      -e "s/HH/${str_hour}/g"     \
      -e "s/mm/${str_minute}/g"  \
      -e "s/SS/${str_second}/g"  < ./input.nml.template > input.nml

./convert_GDC_ncdf 
mv obs_seq.gdc obs_seq.gdc.${str_year}${str_month}${str_day}${str_hour}
let hour=hour+1

done

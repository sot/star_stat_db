#!/bin/bash

SKA=/proj/sot/ska
if [ -d regress ];
then
    rm -r regress
fi
mkdir regress

for obsid in 52789 10539 4639 9600 13894 16612
do
    ${SKA}/share/star_stat_db/update_star_stats.py --obsid ${obsid} \
             --missing-list ${SKA}/data/star_stat_db/ok_missing.csv \
             --server regress/flight.db3 --verbose 2
    ./update_star_stats.py --obsid ${obsid} --server regress/test.db3 --verbose 2
done
for table in trak_stats_data trak_slot_warnings acq_stats_data acq_stats_warnings
do
    ./write_out_table.py regress/test.db3 ${table} regress/test_${table}.dat
    ./write_out_table.py regress/flight.db3 ${table} regress/flight_${table}.dat
done
for table in trak_stats_data trak_slot_warnings acq_stats_data acq_stats_warnings
do
   diff regress/flight_${table}.dat regress/test_${table}.dat >> regress/regress_diffs
done

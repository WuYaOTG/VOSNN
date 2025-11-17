
set Project     hls_cosim.prj
set Solution    solution_2
set Device      "xcu50-fsvh2104-2L-e"
set Flow        "vitis"
set Clock       5

open_project $Project -reset   

set_top Ds_Oc_KNN_search_hw

add_files DSVS_Octree.cpp -cflags -I.
add_files DSVS_Octree_hw.cpp -cflags -I.
add_files DSVS_Octree.h -cflags -I
add_files -tb main.cpp -cflags -I.
add_files -tb timer.h -cflags -I.

open_solution -flow_target $Flow $Solution -reset  
set_part $Device
create_clock -period $Clock -name default

csim_design

puts "successful!!!  the HLS time is [clock format [clock seconds]]"
exec echo "successful!!! HLS start time at [clock format [clock seconds]]" >> "time_record_hls_run.txt"

exit
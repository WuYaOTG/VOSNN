open_project ./build_temp/link/vivado/vpl/prj/prj.xpr
open_run impl_1
report_power -file ./build_temp/link/vivado/vpl/prj/prj.runs/impl_1/full_util_routed_power.rpt
file copy ./build_temp/link/vivado/vpl/prj/prj.runs/impl_1/full_util_routed_power.rpt ./vitis_build.prj/
file copy ./build_temp/link/vivado/vpl/prj/prj.runs/impl_1/full_util_routed.rpt ./vitis_build.prj/
file copy ./build_temp/link//vivado/vpl/prj/prj.runs/impl_1/kernel_util_routed.rpt ./vitis_build.prj/
exit
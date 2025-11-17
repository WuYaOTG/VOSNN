# set the compiling environment
source /tools/Xilinx/Vitis/2024.2/settings64.sh
source /opt/xilinx/xrt/setup.sh
export PLATFORM_REPO_PATHS=/opt/xilinx/platforms/xilinx_u50_gen3x16_xdma_5_202210_1/
export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu

# run the code
time1=$(date)
echo "build host start time at $time1"
make all > build.log 2>&1
tail -n 2 build.log
time1=$(date)
echo "Finish host end time at $time1"

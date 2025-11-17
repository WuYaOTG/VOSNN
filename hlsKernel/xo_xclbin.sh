# set the compiling environment
source /tools/Xilinx/Vivado/2024.2/settings64.sh
source /opt/xilinx/xrt/setup.sh
export PLATFORM_REPO_PATHS=/opt/xilinx/platforms/xilinx_u50_gen3x16_xdma_5_202210_1/
export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu

# run the code
time1=$(date)
echo "build kernel start time at $time1"

make TARGET=hw > build.log 2>&1
tail -n 2 build.log

time1=$(date)
echo "Finish kernel end time at $time1"
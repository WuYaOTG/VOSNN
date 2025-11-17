
## VOSNN

An Efficient Bitonic-Based Voxel-Enhanced Octree Nearest Neighbor Search Accelerator for Large-Scale Point Cloud (for Vitis HLS)

## Overview

VOSNN is a high-performance accelerator design targeting large-scale point cloud nearest neighbor search. It leverages a bitonic sorting based approach combined with voxel-enhanced octree data structures to optimize nearest neighbor search. The project is implemented using Vitis HLS 2024.2 and consists of both PS (Processing System) and PL (Programmable Logic) components.

## Environment Setup

Since VOSNN uses Vitis HLS 2024.2, first source the environment:

```bash
source /tools/Xilinx/Vitis/2024.2/settings64.sh
```

Set unlimited stack size:

```bash
ulimit -s unlimited
```

Make sure you have sufficient memory and system resources.

## Project Structure

- **src/**  
  Contains the core HLS C/C++ source code.

- **hlsKernel/**  
  Scripts for converting generated `.xo` kernel files to `.xclbin` format.

- **host/**  
  Host-side code to interface with the FPGA and run the accelerator.

## Build and Run Instructions

1. Navigate to the core source directory:

```bash
cd src
```

2. Run C simulation (CSIM):

```bash
vitis-run --mode hls --tcl run_csim.tcl
```

3. Run co-simulation (COSIM) or generate `.xo` kernel files:

```bash
vitis-run --mode hls --tcl run_hw.tcl
```

4. Create `.xclbin` file for FPGA deployment:

```bash
cd ../hlsKernel
./xo_xclbin.sh
```

Wait until the `.xclbin` file is generated.

5. Compile and run the host program:

```bash
cd ../host
./run_host.sh
```

This will compile `host.cpp` and launch the host application which loads the `.xclbin` onto the FPGA and runs the accelerator.

## Notes

- The PS-PL communication and accelerator control are handled within the host application.
- Make sure all prerequisites for Vitis HLS and FPGA deployment are installed and configured properly.
- Sufficient system memory and stack size are critical for successful simulation and synthesis.

## Contact

For questions or issues, please contact the project maintainer.
```
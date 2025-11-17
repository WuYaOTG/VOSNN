BUILD_DIR     := ./build_temp
SRCDIR       := /home/wyt/workspack/vitis/DsOc_v2/src

SOURCES := host.cpp DSVS_Octree.cpp DSVS_Octree_hw.cpp
HOST_OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

HOST_EXE     := $(BUILD_DIR)/host.exe

HOST_CXX ?= g++

CXXFLAGS += -I$(SRCDIR)
CXXFLAGS += -I$(XILINX_XRT)/include 
CXXFLAGS += -I$(XILINX_HLS)/include 
CXXFLAGS += -I$(XILINX_VIVADO)/include 
CXXFLAGS += -I/usr/include 
CXXFLAGS += -std=c++17 -O0 -g -Wall -Wno-unknown-pragmas -Wno-unused-label -Wno-narrowing -fmessage-length=0
CXXFLAGS += -DVITIS_PLATFORM=$(VITIS_PLATFORM)
CXXFLAGS += -c

LDFLAGS += -lxrt_coreutil -lpthread
LDFLAGS += -L$(XILINX_XRT)/lib

.PHONY: all clean

all: $(HOST_EXE)

clean:
	rm -rf $(BUILD_DIR) *.log $(HOST_EXE)

$(BUILD_DIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(HOST_CXX) $(CXXFLAGS) -o $@ $<

$(HOST_EXE): $(HOST_OBJECTS)
	@echo "Linking Host executable with objects: $^"
	$(HOST_CXX) -o $@ $^ $(LDFLAGS)
	@echo "GENERATE $@ Done!"


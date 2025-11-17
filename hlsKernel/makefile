# Makefile for Vitis Accelerated Application

# 用户配置区域
# PLATFORM      := xilinx_zcu104_base_202420_1
PLATFORM      := xilinx_u50_gen3x16_xdma_5_202210_1
XO_FILE       := ./Ds_Oc_KNN_search_hw.xo
XCLBIN_NAME   := DsOc_kernel
CONFIG_FILE   := ./link.cfg
TARGET        ?= hw            # 支持 hw/hw_emu/sw_emu
ENABLE_PROFILING ?= 0          # 启用性能分析（1:启用，0:禁用）

# 目录配置
BUILD_DIR     := ./build_temp
REPORT_DIR    := ./report
OUTPUT_DIR    := ./output

# 工具配置
VPP           := v++
XCLBIN        := $(OUTPUT_DIR)/$(XCLBIN_NAME).xclbin

# 根据目标模式设置参数 ==============================================
VPP_FLAGS := -l
VPP_FLAGS += -t $(TARGET)      # 设置构建目标模式
# 性能分析配置
ifeq ($(ENABLE_PROFILING),1)
VPP_FLAGS += --profile_kernel data:all:all:all
VPP_FLAGS += --profile.data all
endif

# 调试配置（硬件仿真时启用波形）
ifeq ($(TARGET),hw_emu)
VPP_FLAGS += --debug --optimize 0
VPP_FLAGS += --profile.data all:all:all
endif

# 默认目标
all: xclbin

# 主构建目标
xclbin: $(XCLBIN)

$(XCLBIN): $(XO_FILE) $(CONFIG_FILE)
# 创建必要目录
	mkdir -p $(BUILD_DIR) $(REPORT_DIR) $(OUTPUT_DIR)
# 执行链接命令
	$(VPP) -l \
		--platform $(PLATFORM) \
		--temp_dir $(BUILD_DIR) \
		--report_dir $(REPORT_DIR) \
		--config $(CONFIG_FILE) \
		--save-temps \
		-o $@ \
		$(XO_FILE) 2>&1 | tee $(BUILD_DIR)/link.log
	@echo "XCLBIN generated: $@"

# 辅助目标
check-env:
	@which $(VPP) || (echo "Vitis environment not configured!"; exit 1)

clean:
	rm -rf $(BUILD_DIR) $(REPORT_DIR) $(OUTPUT_DIR)
	@echo "Cleaned build artifacts"

# 帮助信息
help:
	@echo "Usage:"
	@echo "	make [TARGET=<hw|hw_emu|sw_emu>] [ENABLE_PROFILING=1]"
	@echo ""
	@echo "Examples:"
	@echo "	make TARGET=hw_emu ENABLE_PROFILING=1
	@echo "	make TARGET=sw_emu

.PHONY: all xclbin check-env clean
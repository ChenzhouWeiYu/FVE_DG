# common.mk - 公共配置文件
# 必须放在项目根目录下
# 被各 example 目录下的 Makefile include
# 自动探测自身位置，并设置 ROOT_DIR

# 获取当前 common.mk 的绝对路径（Makefile 中第一个被读取的文件）
ifndef COMMON_MK_ALREADY_INCLUDED
COMMON_MK_ALREADY_INCLUDED := 1

# 获取 common.mk 文件的真实路径
override THIS_MAKEFILE := $(realpath $(lastword $(MAKEFILE_LIST)))
override COMMON_MK := $(THIS_MAKEFILE)
override ROOT_DIR := $(dir $(COMMON_MK))

# 去除末尾斜杠，防止路径拼接错误
override ROOT_DIR := $(patsubst %/,%,$(ROOT_DIR))

# $(info [DEBUG] THIS_MAKEFILE = $(THIS_MAKEFILE))
# $(info [DEBUG] COMMON_MK = $(COMMON_MK))
# $(info [DEBUG] ROOT_DIR = $(ROOT_DIR))
# $(info [DEBUG] INCLUDE_FLAGS = $(INCLUDE_FLAGS))

# 编译器与标志
CXX := g++
CXXFLAGS := -std=c++17 -O3 -march=native
CXXFLAGS += -Wall -Wno-unused-variable -Wno-unused-but-set-variable -Wno-comments
CXXFLAGS += -fopenmp -mfma -mavx2 -DCGAL_HEADER_ONLY -DCGAL_DISABLE_GMP
LDFLAGS := -fopenmp



# 包含路径
INCLUDE_FLAGS := -I$(ROOT_DIR)/include -I$(ROOT_DIR)/external

# 构建输出目录（所有 src/.o 都放这里）
BUILD_DIR := $(ROOT_DIR)/build
OBJ_DIR := $(BUILD_DIR)/obj

# src 路径
SRC_DIR := $(ROOT_DIR)/src

# 自动扫描 src 下所有 .cpp 文件
SRC_SUBDIRS := $(shell find $(SRC_DIR) -type d)
SHARED_SRCS := $(foreach dir,$(SRC_SUBDIRS),$(wildcard $(dir)/*.cpp))

# 对应的 .o 文件路径（build/obj/src/...）
SHARED_OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SHARED_SRCS))

# 工具函数：创建目录
MKDIR_P := mkdir -p

# 共享源文件的编译规则
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@


endif # COMMON_MK_ALREADY_INCLUDED
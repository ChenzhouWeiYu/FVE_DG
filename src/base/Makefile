# src/Mesh/Makefile

# 基础配置
CXX := g++
CXXFLAGS := -std=c++17 -O3 -g 
CXXFLAGS += -fopenmp -mfma -mavx2 -DCGAL_HEADER_ONLY -DCGAL_DISABLE_GMP
ROOT_DIR := $(shell until [ -d "include" -o "$$PWD" = "/" ]; do cd ..; done; echo $$PWD)
INCLUDE_FLAGS := -I$(ROOT_DIR)/include -I$(ROOT_DIR)/external
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.d
DEP_DIR := .deps

# 文件发现
SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)

# 默认目标
all: $(OBJS)

# 编译规则
%.o: %.cpp $(DEP_DIR)/%.d | $(DEP_DIR)
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(INCLUDE_FLAGS) -c $< -o $@

# 依赖目录
$(DEP_DIR):
	@mkdir -p $@

# 包含依赖
$(DEP_DIR)/%.d: ;
.PRECIOUS: $(DEP_DIR)/%.d
-include $(patsubst %,$(DEP_DIR)/%.d,$(basename $(SRCS)))

clean:
	rm -f *.o
	rm -rf $(DEP_DIR)

.PHONY: all clean
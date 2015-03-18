PROGRAM_NAME := program

RM := rm -rf
DATE := date
SRC_PATH := src
BUILD_PATH := build
BIN_PATH := bin
SRC_EXT := cpp
COMPILER := mpiCC
STANDART := -std=c++0x
HOSTFILE := hosts
HOSTS := --hostfile ${HOSTFILE}
L_FLAGS := -lm -lrt 

# Macros for timing compilation
TIME_FILE = $(dir $@).$(notdir $@)_time
START_TIME = $(DATE) '+%s' > $(TIME_FILE)
END_TIME = read st < $(TIME_FILE) ; \
	$(RM) $(TIME_FILE) ; \
	st=$$((`$(DATE) '+%s'` - $$st - 86400)) ; \
	echo `$(DATE) -u -d @$$st '+%H:%M:%S'`

# Verbose option, to output compile and link commands
export V := false
export CMD_PREFIX := @
ifeq ($(V),true)
	CMD_PREFIX :=
endif

WARNINGS_ERRORS := #-pedantic -Wall -Wextra -Wno-deprecated -Wno-unused-parameter  -Wno-enum-compare -Weffc++

debug: COMPILER += -DDEBUG -g
debug: export EXCLUDED_MAIN_FILE := tests.cpp
debug: export EXCLUDED_DIRECTORY := */tests/*
debug: export BUILD_PATH := build/debug
debug: export BIN_PATH := bin/debug

release: COMPILER += -O3
release: export EXCLUDED_MAIN_FILE := tests.cpp
release: export EXCLUDED_DIRECTORY := */tests/*
release: export BUILD_PATH := build/release
release: export BIN_PATH := bin/release

test: export EXCLUDED_MAIN_FILE := main.cpp
test: export EXCLUDED_DIRECTORY :=
test: export BUILD_PATH := build/test
test: export BIN_PATH := bin/test

install: export BIN_PATH := bin/release

SRC_FILES := $(shell find $(SRC_PATH)/ -name '*.$(SRC_EXT)' ! -iname '$(EXCLUDED_MAIN_FILE)' -not -path '$(EXCLUDED_DIRECTORY)' | sort -k 1nr | cut -f2-)
OBJS := $(SRC_FILES:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)
DEP := $(OBJS:.o=.d)

################################################################################

.PHONY: release
release: dirs
	@echo "Beginning release build"
	@$(START_TIME)
	@$(MAKE) all --no-print-directory
	@echo "Total build time: "
	@$(END_TIME)

# Debug build for gdb debugging
.PHONY: debug
debug: dirs
	@echo "Beginning debug build"
	@$(START_TIME)
	@$(MAKE) all --no-print-directory
	@echo "Total build time: "
	@$(END_TIME)

# Test build for gtests
.PHONY: test
test: dirs
	@echo "Beginning test build"
	@$(START_TIME)
	@$(MAKE) all --no-print-directory
	@echo "Total build time: "
	@$(END_TIME)

.PHONY: dirs
dirs:
	@echo "Creating directories"
	@mkdir -p $(dir $(OBJS))
	@mkdir -p $(BIN_PATH)
	@echo "Directories created"

# Removes all build files
.PHONY: clean
clean:
	@echo "Deleting $(PROGRAM_NAME) symlink"
	@$(RM) $(PROGRAM_NAME)
	@echo "Deleting directories"
	@$(RM) -r build
	@$(RM) -r bin

.PHONY: all
all: $(BIN_PATH)/$(PROGRAM_NAME)
	@echo "Making symlink: $(PROGRAM_NAME) -> $<"
	@$(RM) $(PROGRAM_NAME)
	@ln -s $(BIN_PATH)/$(PROGRAM_NAME) $(PROGRAM_NAME)

.PHONY: run
run:
	mpirun $(HOSTS) $(PROGRAM_NAME)

################################################################################

$(BIN_PATH)/$(PROGRAM_NAME): $(OBJS)
	@echo 'Linking target: $@'
	@echo 'Invoking: $(NVCC) Linker'
	$(COMPILER) $(STANDART) $(L_FLAGS) -o $(BIN_PATH)/$(PROGRAM_NAME) $(OBJS)
	chmod +x $(BIN_PATH)/$(PROGRAM_NAME)
	@echo 'Finished building target: $@'
	@echo ' '

$(BUILD_PATH)/%.o: $(SRC_PATH)/%.$(SRC_EXT)
	@$(START_TIME)
	@echo 'Building file: $< -> $@'
	@echo 'Invoking: $(COMPILER) Compiler'
	$(COMPILER) $(STANDART) $(DEFINES) $(WARNINGS_ERRORS) $(FLAGS) -c -MMD -MP -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo '\t Compile time: '
	@$(END_TIME)
	@echo ' '

# Top-level build for the LEAP core binaries.

CONDA_GCC := $(firstword $(wildcard $(CONDA_PREFIX)/bin/*-gcc))
CONDA_GXX := $(firstword $(wildcard $(CONDA_PREFIX)/bin/*-g++))

CC ?= $(if $(CONDA_GCC),$(CONDA_GCC),gcc)
CXX ?= $(if $(CONDA_GXX),$(CONDA_GXX),g++)
LD ?= ld

OPT_FLAGS ?= -O3
STD_FLAGS ?= --std=c++11
ISA_FLAGS ?= -mbmi -mavx2 -msse4.2
WARN_FLAGS ?= -fpermissive
LD_RFLAGS ?= -r
USE_PARASAIL ?= 0

NW_PATH = ./needleman_wunsch-0.3.5
LIBS_PATH = $(NW_PATH)/libs
UTILITY_LIB_PATH := $(LIBS_PATH)/utility_lib
STRING_BUF_PATH := $(LIBS_PATH)/string_buffer
BIOINF_LIB_PATH := $(LIBS_PATH)/bioinf
SCORING_PATH := $(LIBS_PATH)/alignment_scoring

CPPFLAGS += -I . \
	-I $(UTILITY_LIB_PATH) \
	-I $(STRING_BUF_PATH) \
	-I $(BIOINF_LIB_PATH) \
	-I $(SCORING_PATH) \
	-I $(NW_PATH) \
	-DCOMPILE_TIME='"$(shell date)"' \
	-DSCORE_TYPE='int'

ifeq ($(origin CONDA_PREFIX),environment)
CPPFLAGS += -I$(CONDA_PREFIX)/include
LDFLAGS += -L$(CONDA_PREFIX)/lib
endif

CXXFLAGS += $(OPT_FLAGS) $(STD_FLAGS) $(ISA_FLAGS) $(WARN_FLAGS)

CORE_EXECUTABLES = popcount bit_convert vectorED testLV testLV_BAG vectorSHD_ED testRefDB
PARASAIL_EXECUTABLES = testNW
EXECUTABLES = $(CORE_EXECUTABLES)

ifeq ($(USE_PARASAIL),1)
EXECUTABLES += $(PARASAIL_EXECUTABLES)
PARASAIL_LDLIBS = -lparasail
endif

all: $(EXECUTABLES)

smoke: testLV testLV_BAG vectorSHD_ED testRefDB vectorED

LV_BAG.o: LV_BAG.cc LV_BAG.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

LV.o: LV.cc LV.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

SIMD_ED.o: SIMD_ED.cc SIMD_ED.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

print.o: print.c print.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

popcount.o: popcount.c popcount.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

popcount: popcount.o print.o popcountMain.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

bit_convert.o: bit_convert.c bit_convert.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

bit_convert: print.o bit_convert.o bit_convertMain.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

shift.o: shift.c shift.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

SHD.o: SHD.cc SHD.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

mask.o: mask.c mask.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

sse.o: mask.o print.o bit_convert.o popcount.o vector_filter.o
	$(LD) $(LD_RFLAGS) $^ -o $@

vector_filter: mask.o print.o bit_convert.o popcount.o vector_filter.o vector_filterMain.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

countPassFilter: mask.o print.o bit_convert.o popcount.o vector_filter.o countPassFilter.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

timeSSE: timeSSE.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

needleman_wunsch.o: $(NW_PATH)/needleman_wunsch.c $(NW_PATH)/needleman_wunsch.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

vectorED: SIMD_ED.o print.o bit_convert.o vectorED.cc shift.o SHD.o mask.o popcount.o needleman_wunsch.o $(wildcard $(SCORING_PATH)/*.c) $(UTILITY_LIB_PATH)/utility_lib.c $(BIOINF_LIB_PATH)/bioinf.c $(STRING_BUF_PATH)/string_buffer.c $(NW_PATH)/nw_cmdline.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -lz

testNW: SIMD_ED.o print.o bit_convert.o testNW.cc shift.o SHD.o mask.o popcount.o needleman_wunsch.o $(wildcard $(SCORING_PATH)/*.c) $(UTILITY_LIB_PATH)/utility_lib.c $(BIOINF_LIB_PATH)/bioinf.c $(STRING_BUF_PATH)/string_buffer.c $(NW_PATH)/nw_cmdline.c
ifeq ($(USE_PARASAIL),1)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) -lz $(PARASAIL_LDLIBS)
else
	@echo "testNW requires USE_PARASAIL=1 and libparasail." >&2
	@false
endif

testLV: LV.o testLV.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

testLV_BAG: LV_BAG.o testLV_BAG.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

vectorSHD_ED: SIMD_ED.o SHD.o mask.o print.o bit_convert.o shift.o popcount.o vectorSHD_ED.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

test_SIMD_ED: SIMD_ED.o vector_filter.o bit_convert.o mask.o popcount.o print.o test_ED.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

diffED: LV.o SIMD_ED.o mask.o print.o bit_convert.o popcount.o vector_filter.o diffED.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

RefDB.o: RefDB.cc RefDB.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

testRefDB: RefDB.o bit_convert.o shift.o print.o RefDBMain.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@

.PHONY: all clean smoke

clean:
	rm -f $(CORE_EXECUTABLES) $(PARASAIL_EXECUTABLES) *.o

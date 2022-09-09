FC := h5pfc

ifeq (${GCC},y)
# gcc:
FFLAGS := -Ofast -O3 -DSoA -ffree-line-length-512
else ifeq (${INTEL},y)
# intel:
FFLAGS := -Ofast -O3 -DSoA -ipo -qopenmp-simd -qopt-report=5 -qopt-streaming-stores auto
else ifeq (${INTELAVX512},y)
# intel (avx512):
FFLAGS := -Ofast -O3 -DSoA -ipo -qopenmp-simd -qopt-report=5 -qopt-streaming-stores auto -xCORE-AVX512 -qopt-zmm-usage:high
else ifeq (${INTELAVX2},y)
# intel (avx2):
FFLAGS := -Ofast -O3 -DSoA -ipo -qopenmp-simd -qopt-report=5 -qopt-streaming-stores auto -xCORE-AVX2
else ifeq (${AMD},y)
# amd:
FFLAGS := -Ofast -O3 -DSoA -fwhole-program -mavx2 -fopt-info-vec -fopt-info-vec-missed -ftree-vectorizer-verbose=5
else
$(error "unknown compiler")
endif

.PHONY: compile clean

compile : compile_test.o
	${FC} $(FFLAGS) $< -o $@

%.o : %.F90
	${FC} $(FFLAGS) -c $< -o $@

clean :
	rm *.o compile
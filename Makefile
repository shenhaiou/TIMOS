CXX    = clang++
TARGET = source/timos

SRCS = source/mesh.cpp \
       source/optics.cpp \
       source/source_io.cpp \
       source/simulation.cpp \
       source/output.cpp \
       source/main.cpp

# Base flags for M1 Max
CXXFLAGS  = -O3 -std=c++23 -march=native -ffast-math -funroll-loops
CXXFLAGS += -DHAS_ACCELERATE -framework Accelerate
CXXFLAGS += -Isource

# PGO artefacts
PROF_RAW  = source/timos.profraw
PROF_DATA = source/timos.profdata
PROF_EXE  = source/timos_pgo

PROF_ARGS = -f example/onelayer/one_layer_18_18_1_2.mesh \
            -s example/onelayer/one_layer_18_18_1_2.source \
            -p example/onelayer/mua005_mus05.opt \
            -m is -o /tmp/timos_profile_out.dat -t 10

TESTFLAGS = -O1 -std=c++23 -Isource -DHAS_ACCELERATE -framework Accelerate
TESTS     = tests/test_mesh tests/test_optics tests/test_integration
# All modules except main.cpp (each test provides its own main)
LIB_SRCS  = source/mesh.cpp source/optics.cpp source/source_io.cpp \
             source/simulation.cpp source/output.cpp

.PHONY: all release pgo clean test

all: release

release: $(SRCS)
	$(CXX) $(CXXFLAGS) -flto -o $(TARGET) $(SRCS)

pgo: $(SRCS)
	@echo "=== Step 1: build instrumented binary ==="
	$(CXX) $(CXXFLAGS) -fprofile-instr-generate -o $(PROF_EXE) $(SRCS)
	@echo "=== Step 2: collect profile ==="
	LLVM_PROFILE_FILE=$(PROF_RAW) $(PROF_EXE) $(PROF_ARGS)
	@echo "=== Step 3: merge profile ==="
	xcrun llvm-profdata merge -output=$(PROF_DATA) $(PROF_RAW)
	@echo "=== Step 4: build optimised binary ==="
	$(CXX) $(CXXFLAGS) -flto -fprofile-instr-use=$(PROF_DATA) -o $(TARGET) $(SRCS)
	@echo "=== Done: $(TARGET) ==="

test: $(TESTS)
	@echo "=== Running tests (from project root) ==="
	./tests/test_mesh
	./tests/test_optics
	./tests/test_integration
	@echo "=== All tests passed ==="

tests/test_mesh: tests/test_mesh.cpp $(LIB_SRCS)
	$(CXX) $(TESTFLAGS) -o $@ $^

tests/test_optics: tests/test_optics.cpp $(LIB_SRCS)
	$(CXX) $(TESTFLAGS) -o $@ $^

tests/test_integration: tests/test_integration.cpp $(LIB_SRCS)
	$(CXX) $(TESTFLAGS) -o $@ $^

clean:
	rm -f $(TARGET) $(PROF_EXE) $(PROF_RAW) $(PROF_DATA) $(TESTS)

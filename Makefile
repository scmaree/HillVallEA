CC = g++
CFLAGS = -O3 -std=c++11 -MMD

CEC_DIR := ./CEC2013_niching_benchmark
CEC_SRC_FILES := $(wildcard $(CEC_DIR)/*.cpp)
CEC_OBJ_FILES := $(patsubst $(CEC_DIR)/%.cpp,$(CEC_DIR)/%.o,$(CEC_SRC_FILES))
CEC_DEP_FILES := $(patsubst $(CEC_DIR)/%.cpp,$(CEC_DIR)/%.d,$(CEC_SRC_FILES))

HVEA_DIR := ./HillVallEA
HVEA_SRC_FILES := $(wildcard $(HVEA_DIR)/*.cpp)
HVEA_OBJ_FILES := $(patsubst $(HVEA_DIR)/%.cpp,$(HVEA_DIR)/%.o,$(HVEA_SRC_FILES))
HVEA_DEP_FILES := $(patsubst $(HVEA_DIR)/%.cpp,$(HVEA_DIR)/%.d,$(HVEA_SRC_FILES))

all: example_cec2013_benchmark example_simple

example_cec2013_benchmark: example_CEC2013_benchmark.o $(HVEA_OBJ_FILES) $(CEC_OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ example_CEC2013_benchmark.o $(HVEA_OBJ_FILES) $(CEC_OBJ_FILES)

example_simple: example_simple.o $(HVEA_OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ example_simple.o $(HVEA_OBJ_FILES)
	
%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(CEC_OBJ_FILES) $(CEC_DEP_FILES) $(HVEA_OBJ_FILES) $(HVEA_DEP_FILES) *.d *.o

clean_run:
	rm -f example_cec2013_benchmark example_simple elites.dat statistics.dat

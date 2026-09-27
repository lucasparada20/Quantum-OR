CXX = g++

CPLEX_DIR = /some/path/to/Cplex

CXXFLAGS = -O3 \
    -I$(CPLEX_DIR)/cplex/include \
    -I$(CPLEX_DIR)/concert/include

LDFLAGS = \
    -L$(CPLEX_DIR)/cplex/lib/x86-64_linux/static_pic \
    -L$(CPLEX_DIR)/concert/lib/x86-64_linux/static_pic

LDLIBS = -lilocplex -lcplex -lconcert -lm -lpthread



TARGET = vrp_qubo
SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) $(LDFLAGS) $(LDLIBS) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o


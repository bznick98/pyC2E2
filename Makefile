#Directories
BIN_DIR  = ./bin
EXEC_DIR = ./c2e2/cpp_backend/lib
OBJ_DIR  = ./obj
SRC_DIR  = ./c2e2/cpp_backend/cpp

#Compiler
CXX = g++
CFLAGS = -g -std=c++11
INCLUDES = -I $(SRC_DIR) -I /usr/local/Cellar/python@3.8/3.8.4/Frameworks/Python.framework/Versions/3.8/include/python3.8 -I /usr/local/Cellar/eigen
LFLAGS = -lboost_python38 -lpython3.8 -ldl -lppl -lglpk -lgmp
LIBLDIR = /usr/local/Cellar/python@3.8/3.8.4/Frameworks/Python.framework/Versions/3.8/lib/python3.8/config-3.8-darwin/

CPP_SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
CPP_OBJECTS := $(CPP_SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

#Target
EXEC     = libc2e2.so
TARGET  :=  $(EXEC_DIR)/$(EXEC)

.PHONY: all
all: $(TARGET)

.PHONY: dirs
dirs:
	@ mkdir -p $(OBJ_DIR)
	@ mkdir -p $(BIN_DIR)

# FIXME
$(TARGET): $(CPP_OBJECTS)
	# -Wl,--no-undefined cannot work on MacOS
	$(CXX) -shared -Wl,-undefined,error $(CFLAGS) -o $@ $(CPP_OBJECTS) $(LFLAGS) $(INCLUDES) -L$(LIBLDIR)

# FIXME
$(CPP_OBJECTS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CFLAGS) -fPIC -c -o $@ $< $(INCLUDES)
	#$(CXX) $(CFLAGS) $(INCLUDES) -fPIC -c -o $@ $< 

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR)/*.o
	rm -rf $(EXEC_DIR)/*.so
CXX = g++
WRN = -Weffc++ -Wall -Wextra -Wconversion -Wsign-conversion
STD = -std=c++17
DBG = -g
OPT = -O0
LIB = # -lglut -lGL
CXXFLAGS = $(WRN) $(OPT) $(DBG) $(STD)

SRC_EXT = cpp

SRC_DIR = source
HED_DIR = include
OBJ_DIR = .object

SRC = $(wildcard $(SRC_DIR)/*.$(SRC_EXT))
HED = $(wildcard HED_DIR/*)
OBJ = $(patsubst $(SRC_DIR)%, $(OBJ_DIR)%, $(patsubst %.$(SRC_EXT), %.o, $(SRC)))
OUT = binary

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.$(SRC_EXT) $(HED)
	$(CXX) -c $< $(CXXFLAGS) -o $@ -I$(HED_DIR)

$(OUT) : $(OBJ) 
	$(CXX) $(OBJ) -o $@ $(LIB)

initialize:
	mkdir $(OBJ_DIR)

clean:
	rm $(OBJ_DIR) -r
	rm $(OUT)


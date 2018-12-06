# Makefile
# Author     : Yueze Tan
# Email      : yut75@psu.edu
# Written for: CSE 597-002(Fall 2018), HW01
# Last modified @ 2018/09/22

NAME = LUDecomposition

CC = g++
CFLAGS = -O3 -mfma -mavx -std=c++11 -Wno-unused-result
LFLAGS = -lm

SRCPATH = ./src
OBJPATH = ./build
BINPATH = ./bin

SOURCES = $(wildcard ${SRCPATH}/*.cpp)
OBJECTS = $(patsubst %.cpp, ${OBJPATH}/%.o, $(notdir ${SOURCES}))
BINNAME = $(NAME)

program: $(OBJECTS) | $(BINPATH)
	$(CC) $(OBJECTS) -o $(BINPATH)/$(BINNAME) $(LFLAGS)
	@echo ""
	@echo "Make completed."
	@echo ""

$(OBJPATH)/%.o: $(SRCPATH)/%.cpp | $(OBJPATH)
	$(CC) -c $< -o $@ $(CFLAGS)

$(BINPATH):
	@mkdir $(BINPATH)

$(OBJPATH):
	@mkdir $(OBJPATH)

clean:
	@rm -rf $(OBJPATH)
	@echo "Clean completed."
	@echo ""


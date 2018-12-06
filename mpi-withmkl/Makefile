# Makefile
# Author     : Yueze Tan
# Email      : yut75@psu.edu
# Written for: CSE 597-002(Fall 2018), HW02
# Last modified @ 2018/10/26

NAME = MKLJacobi

CC = tau_cxx.sh    # for profiling
# CC = mpicxx    # for regular executing
MKL_INCLUDE_PATH = /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/include
MKL_LIB_PATH = /opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/

CFLAGS = -O2 -mfma -mavx -std=c++11 -Wno-unused-result -g -I$(MKL_INCLUDE_PATH) -DMKL_ILP64 -m64
LFLAGS = -lm -Wl,--no-as-needed -L$(MKL_LIB_PATH) -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread

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

# part of this Make file is taken from 
# the sample at lecure notes
FOO := toolsProject.out
SRC := $(wildcard *.cpp) 
OBJ := $(patsubst %.cpp,%.o,$(SRC))
# Options
CC := g++
CFLAGS := -I$(TACC_GSL_INC) -I$(TACC_GRVY_INC) 
LDFLAGS := -L$(TACC_GSL_LIB) -L$(TACC_GRVY_LIB)
LDLIBS := -lgsl -lgslcblas -lm -lgrvy 
# # Rules
$(FOO): $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.out: %.cpp
	$(CC) $(CFLAGS) -c $<
toolsProject.o := header.h
# Useful phony targets
.PHONY: clean neat 
clean: neat 
	$(RM) $(EXEC)
neat: 
	$(RM) $(OBJ)

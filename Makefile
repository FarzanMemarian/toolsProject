# A big part of this Make file is taken from 
# the sample at lecure notes
FOO := main.out
SRC := $(wildcard *.cpp) 
OBJ := $(patsubst %.cpp,%.o,$(SRC))
# Options
CC := g++
#CFLAGS := -I$(TACC_GSL_INC) # from build_me file
#LDFLAGS := -L$(TACC_GSL_LIB) # from build_me file
LDLIBS := -lgsl -lcblas -lm # from build_me file
# # Rules
$(FOO): $(OBJ)
	$(CC) $^ $(LDLIBS) -o $@
%.out: %.cpp
	$(CC) -c $<
main.o := header.h
# Useful phony targets
.PHONY: clean neat 
clean: neat 
	$(RM) $(EXEC)
neat: 
	$(RM) $(OBJ)

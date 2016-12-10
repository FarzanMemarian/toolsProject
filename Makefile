#Files
FOO := toolsProject
SRC := $(wildcard *.cpp)
OBJ := $(patsubst %.cpp,%.o,$(SRC))
# Options
CC := g++
CFLAGS := -I$(TACC_GSL_INC) -I$(TACC_GRVY_INC)
LDFLAGS := -L$(TACC_GSL_LIB) -L$(TACC_GRVY_LIB)
LDLIBS := -lgrvy -lgsl -lgslcblas -lm
# Rules
$(FOO): $(OBJ)
	$(CC) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<
toolsProject := headers.h
check1:
	./ input1.dat > output1.dat
	diff knownAns1.dat output1.dat
#test2:
#	./ODE input2.txt > output2.txt
#	diff correct2.txt output2.txt
check: check1  
# Useful phony targets
.PHONY: clean neat
clean: neat
	$(RM) $(EXEC)
neat:
	$(RM) $(OBJ)

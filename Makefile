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
toolsProject := header.h
#test1:
#	./ODE input1.txt > output1.txt
#	diff knownAns.txt output1.txt
#test2:
#	./ODE input2.txt > output2.txt
#	diff correct2.txt output2.txt
#check: test1 test2 
# Useful phony targets
.PHONY: clean neat
clean: neat
	$(RM) $(EXEC)
neat:
	$(RM) $(OBJ)

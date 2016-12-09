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
#	diff correct1.txt output1.txt
#test2:
#	./ODE input2.txt > output2.txt
#	diff correct2.txt output2.txt
#test3:
#	./ODE input3.txt > output3.txt
#	diff correct3.txt output3.txt
#test4:
#	./ODE input4.txt > output4.txt
#	diff correct4.txt output4.txt
#test5:
#	./ODE input5.txt > output5.txt
#	diff correct5.txt output5.txt
#check: test1 test2 test3 test4 test5
# Useful phony targets
.PHONY: clean neat
clean: neat
	$(RM) $(EXEC)
neat:
	$(RM) $(OBJ)

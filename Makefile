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
# the output of the myEuler function is compared here
test1: 
	./toolsProject input1.dat 
	diff prob1_MyEuler_KNOWN.dat prob1_MyEuler.dat
# the output of the gsl for problem 1 is compared here
test2: 
	./toolsProject input2.dat 
	diff prob1_GSLSolver_KNOWN.dat prob1_GSLSolver.dat 
# the output of the gsl rk2 is compared here
test3: 
	./toolsProject input3.dat
	diff prob2_rk2_KNOWN.dat prob2_rk2.dat 
# the output of the gsl rk4 is compared here
test4: 
	./toolsProject input4.dat
	diff prob2_rk4_KNOWN.dat prob2_rk4.dat 
# the output of the gsl rkf45 is compared here
test5: 
	./toolsProject input5.dat
	diff prob2_rkf45_KNOWN.dat prob2_rkf45.dat 
check: test1 test2 test3 test4 test5 
# Useful phony targets
.PHONY: clean neat
clean: neat
	$(RM) $(EXEC)
neat:
	$(RM) $(OBJ)

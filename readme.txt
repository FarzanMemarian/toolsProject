to run the code on stampede, first run
module load gcc/4.7.1
module load gsl
module load grvy

to run the code, go to the src folder and run

make
make clean
make

this will create the executable named toolsProject

then type:

./toolsProject input.dat

input.dat is an input file which is in the same folder. the choice between the problems and other choices can be made in this file. 

if you want to perform make check, just type:
make check

this will read the input files from the test directory and compare the output of running them with known output

to plot the results, you can run the python file. in a folder where data are available, such as the error folder. note that the outputs do not go directly into error folder and need to be done manually. 

CC = g++

FLAGS = -O3 -std=c++11 -fopenmp

dirLib = $(PWD)/../../../
dirH = $(dirLib)Headers/Lib/
dirHO = $(dirLib)Headers/Objects/
dirNum = $(PWD)/../../../../Numerical/
dirEH = $(dirLib)External/EH/

varDIR = -DDIRCOSMO=\"$(dirLib)\" -DDIRL=\"$(PWD)/\"

FLAGS_LIB = -Wl,-rpath,$(dirLib) -L$(dirLib) -lFUNC -lSTAT -lCOSM -lCM -lCAT -lLN -lRANDOM -lTWOP -lTHREEP -lGLOB
FLAGS_INC = -I$(dirH) -I$(dirHO) -I$(dirNum) -I$(dirEH) 

OBJ = chi2.o
OBJ1 = likelihood.o

ES = so

ifeq ($(SYS),MAC)
        ES = dylib
endif

chi2: $(OBJ) 
	$(CC) $(OBJ) -o chi2 $(FLAGS_LIB) $(FLAGS_GSL)

likelihood: $(OBJ1) 
	$(CC) $(OBJ1) -o likelihood $(FLAGS_LIB) $(FLAGS_GSL)

clean:
	rm *.o chi2 *~ \#* temp* core* -f

chi2.o: chi2.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c chi2.cpp

likelihood.o: likelihood.cpp makefile $(dirLib)*.$(ES) 
	$(CC) $(FLAGS) $(FLAGS_INC) $(varDIR) -c likelihood.cpp

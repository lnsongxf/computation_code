% Compile mex file
mex hello_mex.cpp -I"../CRoutines"

% Openmp version
% mex -v hello_mex.cpp -I"../CRoutines" COMPFLAGS="$COMPFLAGS /Openmp" -DOMP
% Compile mex file
mex vfi_mex.cpp -I"../CRoutines" -DOMP -DBZ_THREADSAFE COMPFLAGS="$COMPFLAGS /Qopenmp" OPTIMFLAGS="/Z7 /O3"
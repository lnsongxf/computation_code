function COMPILE
MKL_ROOT = 'C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\mkl';
MKL_LIB_FOLDER = [MKL_ROOT '\lib\intel64\'];
MKL_LIB1 = ['"' MKL_LIB_FOLDER 'mkl_intel_lp64.lib"'];
MKL_LIB2 = ['"' MKL_LIB_FOLDER 'mkl_core.lib"'];
MKL_LIB3 = ['"' MKL_LIB_FOLDER 'mkl_intel_thread.lib"'];
MKL_LIB = [' ' MKL_LIB1 ' ' MKL_LIB2 ' ' MKL_LIB3 ' '];
% D:\blitz-0.10\src\globals.cpp 
% COMPILE_STR = ['mex qqsMex.cpp' MKL_LIB '-DBZ_THREADSAFE -DUSE_OMP -I"D:\blitz-0.10" -I"D:\cpp_local\routines" OPTIMFLAGS="/Z7 /O2 /QxAVX /DNDEBUG" COMPFLAGS="$COMPFLAGS /Z7 /Qopenmp -Qopt-report:4" LINKOPTIMFLAGS="/DEBUG"']
COMPILE_STR = ['mex -v HetNKMex.cpp' ' ' '-DBZ_THREADSAFE -DUSE_OMP -I"C:\Users\Wenlan\Dropbox\CRoutines" -I"C:\Users\Wenlan\Dropbox\CRoutines" OPTIMFLAGS="/Z7 /O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /Z7 /Qopenmp /Qopt-report-filter:"qqsMex.cpp" /Qopt-report-phase=vec /Qopt-report=5 /Qopt-report-file=stdout" LINKOPTIMFLAGS="/DEBUG"']
% COMPILE_STR = ['mex -v HetNKMex.cpp' ' ' '-DBZ_THREADSAFE -I"C:\Users\Wenlan\Dropbox\CRoutines" -I"C:\Users\Wenlan\Dropbox\CRoutines" OPTIMFLAGS="/Z7 /O2 /DNDEBUG" COMPFLAGS="$COMPFLAGS /Z7 /Qopenmp /Qopt-report-filter:"qqsMex.cpp" /Qopt-report-phase=vec /Qopt-report=5 /Qopt-report-file=stdout" LINKOPTIMFLAGS="/DEBUG"']

COMPILE_DEBUG_STR = ['mex -v -g qqsMex.cpp "C:\Users\Wenlan\Dropbox\CRoutines\src\globals.cpp"' ' ' '-DBZ_THREADSAFE -DMMDEBUG -I"C:\Users\Wenlan Luo\Dropbox\CRoutines" -I"C:\Users\Wenlan\Dropbox\CRoutines" OPTIMFLAGS="/Z7 /O2 /QxAVX /DNDEBUG" COMPFLAGS="$COMPFLAGS /Z7 /Qopenmp /Qopt-report-filter:"qqsMex.cpp" /Qopt-report-phase=vec /Qopt-report=5 /Qopt-report-file=stdout" LINKOPTIMFLAGS="/DEBUG"']


eval(COMPILE_STR);
% eval(COMPILE_DEBUG_STR);

% mex -g qqsMex.cpp -DBZ_DEBUG D:\blitz-0.10\src\globals.cpp -I"D:\blitz-0.10" -I"D:\cpp_local\routines" -I"C:\Program Files (x86)\Intel\Composer XE\mkl\include"
% -DBZ_DEBUG 
end

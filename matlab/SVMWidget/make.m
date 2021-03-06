% To be used if cmake from the root directory fails OR on windows.

disp('Building SVMWidget object files...');

mex -c -largeArrayDims -O  -I. -I../../include mxUtil.cpp ../../src/asvm.cpp ../../src/asvmdata.cpp ...
     ../../src/util.cpp ../../src/svm.cpp

    disp('Building mxSMOSolver...');
    try
if(isunix)
mex   -largeArrayDims -O  -I. -I../../include mxSMOSolver.cpp ../../src/asvm_smo_solver.cpp asvm.o asvmdata.o util.o svm.o mxUtil.o
elseif(ispc)
mex   -largeArrayDims -O  -I. -I../../include mxSMOSolver.cpp ../../src/asvm_smo_solver.cpp asvm.obj asvmdata.obj util.obj svm.obj mxUtil.obj
end
    catch ex
        warning('SMO solver mexFunction could not be built. Make sure the main library is built successfully');
    end


if(true)
    disp('Building mxNLOPTSolver...');
    try
if(isunix)
            mex  -largeArrayDims -L/usr/local/lib -lnlopt -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I/usr/local/include -I. -I../../include mxNLOPTSolver.cpp ../../src/asvm_nlopt_solver.cpp asvm.o asvmdata.o svm.o util.o mxUtil.o
elseif(ispc)
 mex  -largeArrayDims -L/usr/local/lib -lnlopt -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I/usr/local/include -I. -I../../include mxNLOPTSolver.cpp ../../src/asvm_nlopt_solver.cpp asvm.obj asvmdata.obj svm.obj util.obj mxUtil.obj
end
    catch ex
        warning('NLOPT solver mexFunction could not be built. Make sure USE_NLOPT is set TRUE while building the main library');
    end
    
      disp('Building mx_quadprog solver...');
    try
if(isunix)
            mex  -largeArrayDims -L/usr/local/lib -lnlopt -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I/usr/local/include -I. -I../../include mx_quadprog.cpp ../../src/quadprog_nlopt_solver.cpp
elseif(ispc)
 mex  -largeArrayDims -L/usr/local/lib -lnlopt -O  CFLAGS='-Wall' LDOPTIMFLAGS='-Wl,--rpath -Wl,/usr/local/lib'...
 -I/usr/local/include -I. -I../../include mx_quadprog.cpp ../../src/quadprog_nlopt_solver.cpp 
end
    catch ex
        warning('quadprog solver mexFunction could not be built.');
    end
    
end

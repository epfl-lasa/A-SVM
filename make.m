
disp('Compiling mex files');
cd matlab/libsvm
try 
    make
catch ex
   warning('Cannot compile libsvm mex files!!'); 
end
cd ../..
   
cd matlab/SVMWidget
try 
    make
catch ex
   warning('Cannot compile SVMWidget mex files!!'); 
end

cd ../SVMUtil
try 
    make
catch ex
   warning('Cannot compile SVMWidget mex files!!'); 
end
cd ../..

disp('Done.');
disp('Run <svm_widget> to start application');






Augmented-SVM source code: version 1.0 
=======================
Issued - Nov 1, 2012   
Updated - Nov 2, 2016 (plotting svm doesn't work in R2015a -- need to fix)

This package contains the algorithm for learning the Augmented-SVM 
classifier function for combining multiple non-linear dynamics. 
The algorithm was presented in the paper:

Shukla, A. and Billard, A. ["Augmented-SVM: Automatic space partitioning 
for combining multiple non-linear dynamics."](https://papers.nips.cc/paper/4654-augmented-svm-automatic-space-partitioning-for-combining-multiple-non-linear-dynamics.pdf) Neural Information 
Processing Systems (NIPS) 2012. Tahoe, Nevada.


<p align="center">
<img src="https://github.com/epfl-lasa/A-SVM/blob/master/img/class_1.png" width="390"><img src="https://github.com/epfl-lasa/A-SVM/blob/master/img/class_2.png" width="390">
</p>

##Package Structure
This package is organised as follows. 
- The root folder contains cpp source  files in ```src/``` and ```include/``` which are compiled using cmake into a library
``lib/ASVMLearning.so`` This is then linked to the executable ```bin/train``` which
is a command-line interface to the learning algorithm. 
- The ```matlab/``` folder contains mex interface to the main library and the function ```svmtrain``` from ```libsvm```.
- The ASVM utility functions are in ```matlab/SVMUtil```.
- A GUI is provided in ```matlab/SVMWidget``` to draw trajectories using a mouse or stylus and call the various learning algorithms to get the A-SVM model. This is the calling
point for all the ```cpp/mex``` functions compiled before. 


##Installation
###LINUX
Extract the ASVMLearning folder to any location.
```
>> cd <A-SVM_root_dir>
>> mkdir build
>> cd build
>> ccmake ..
>> make
```

You may optionally choose NLOPT at the ccmake gui. You must have nlopt "make install"-ed
somewhere on you filesystem. If you chose standard install location for NLOPT, ccmake 
will find the required files automatically. If you installed NLOPT in some other location
you may need to manually give the location of the nlopt library and nlopt.hpp in the ccmake 
gui. Also, in any case, do not forget to add the folder containing the nlopt library to 
LD_LIBRARY_PATH.

You can also choose to switch off the NLOPT option and in that case you will not be able to
use the NLOPT algorithms from the matlab gui.

Once ASVMLearning.so is built, you need to compile the mex files from matlab. For that just run
the following command in your matlab command line
```
>> setup_path.m 
>> make 
```
This will take care of compiling all the mex interfaces and adding relevant folders to matlab path. If you get GLIBC_XXX errors on calling the mex functions, you may try to run matlab using the provided script "run_matlab.sh".

To run the matlab widget call:
```
>> svm_widget.m 
```

###WINDOWS
Extract the ASVMLearning folder to any location.
```
>> Create a build directory inside ASVMLearning
>> run cmake-gui with source_folder=ASVMLearning and build_folder=ASVMLearning/build
>> configure
>> Optionally choose NLOPT (See above LINUX installation for details).
>> generate
```
You must compile in the "Release" version otherwise the mex linkage will fail later on.

Once ASVMLearning.so is built, you need to compile the mex files from matlab. For that just run
matlab from the ASVMLearning root folder and run "setup_path.m". This will take care of compiling
all the mex interfaces and adding relevant folders to matlab path.


##ThirdParty
###Libsvm
This package uses the function svmtrain from Libsvm. It contains a modified
subset of the original libsvm source which additionally returns the indices 
of the chosen support vectors within the model. Full version of Libsvm is 
available at 
[http://www.csie.ntu.edu.tw/~cjlin/libsvm](http://www.csie.ntu.edu.tw/~cjlin/libsvm)
Please read the COPYRIGHT file before using Libsvm.

###NLOPT
This package gives access to using all the NLOPT algorithms for learning the
A-SVM model. NLOPT is available at
http://ab-initio.mit.edu/wiki/index.php/NLopt

###IPOPT
This package can optionally use the IPOPT solver if it is compiled with its
MATLAB interface. IPOPT is available at
https://projects.coin-or.org/Ipopt


Current Maintainer: [Nadia Figueroa](http://lasa.epfl.ch/people/member.php?SCIPER=238387)



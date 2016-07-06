
echo ... please set matlab root in this script before it runs
MATLAB_ROOT=/Applications/MATLAB_R2016b.app/

echo .. run with sudo if you haven't done that

echo cd'ing into intel tools folder
cd /opt/intel/mkl/tools/builder
cat blas_example_list > blas_lapack_list
cat lapack_example_list >> blas_lapack_list

echo making libraries
make libuni interface=ilp64 export=blas_lapack_list name=libsingle_mkl_ilp64 threading=sequential
cp libsingle_mkl_ilp64* $MATLAB_ROOT/extern/lib/maci64

cp ../Dependencies/mmx_package/*.{c,cpp,h} $MATLAB_ROOT/extern/lib/maci64/

echo SUCCESS!


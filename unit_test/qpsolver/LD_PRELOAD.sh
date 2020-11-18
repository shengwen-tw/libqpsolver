MKL_LIB_PATH=/opt/intel/mkl/lib/intel64/

MKL_LIBS=$MKL_LIB_PATH/libmkl_def.so
MKL_LIBS=$MKL_LIBS:$MKL_LIB_PATH/libmkl_avx2.so
MKL_LIBS=$MKL_LIBS:$MKL_LIB_PATH/libmkl_core.so
MKL_LIBS=$MKL_LIBS:$MKL_LIB_PATH/libmkl_intel_lp64.so
MKL_LIBS=$MKL_LIBS:$MKL_LIB_PATH/libmkl_sequential.so

export LD_PRELOAD=$MKL_LIBS

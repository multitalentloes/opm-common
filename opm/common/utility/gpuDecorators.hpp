/*
  Copyright 2024 SINTEF Digital

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
    This file will define the macro OPM_HOST_DEVICE to be empty if
    we do not compile for GPU architectures, otherwise we will
    set it to "__device__ __host__" to decorate functions that can
    be called from a hip/cuda kernel
*/

// HAVE_CUDA and USE_HIP are found in config.h
#ifndef OPM_GPUDECORATORS_HPP
  #define OPM_GPUDECORATORS_HPP
  #include <config.h>

  #define STRINGIZE(x) #x
  #define PRINT_MACRO_VALUE(x) #x " = " STRINGIZE(x)

  #if HAVE_CUDA // if we will compile with GPU support

    #if USE_HIP // if we compile for AMD architectures
      #include <hip/hip_runtime.h>
    #else // if we compile for Nvidia architectures
      #include <cuda_runtime.h>
    #endif // END USE_HIP

    // define host and device function attributes
    #define OPM_HOST_DEVICE __device__ __host__
    #define OPM_DEVICE __device__
    #define OPM_HOST __host__
    #define OPM_USING_GPU LITERALLY ANYTHING ELSE
    #pragma message(PRINT_MACRO_VALUE(OPM_USING_GPU))
    #pragma message(PRINT_MACRO_VALUE(OPM_HOST_DEVICE))

    // Define OPM_DEVICE_IF_GPUCC based on whether we are using a GPU compiler
    #if defined(__NVCC__) | defined(__HIPCC__)
      #define OPM_DEVICE_IF_GPUCC __device__
    #else
      #define OPM_DEVICE_IF_GPUCC
    #endif // END IF GPUCC

  #else // if we are not using CUDA/HIP, let the macro be empty
    #define OPM_HOST_DEVICE
    #define OPM_DEVICE
    #define OPM_DEVICE_IF_GPUCC
    #define OPM_HOST
    #undef OPM_USING_GPU
  #endif // END ELSE

#endif // END HEADER GUARD

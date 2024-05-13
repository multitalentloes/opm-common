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
#include <config.h>

#if HAVE_CUDA // if we will compile with GPU support

#if USE_HIP // if we compile for AMD architectures
#include <hip/hip_runtime.h>
#else // if we compile for Nvidia architectures
#include <cuda_runtime.h>
#endif

// define host and device function attributes
#define OPM_HOST_DEVICE __device__ __host__

#else // if we are not using CUDA/HIP, let the macro be empty
#define OPM_HOST_DEVICE
#endif

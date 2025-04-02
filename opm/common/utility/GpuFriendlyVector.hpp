/*
  Copyright 2025 EQUINOR

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef OPM_COMMON_UTILITY_GPUFRIENDLYVECTOR_HPP
#define OPM_COMMON_UTILITY_GPUFRIENDLYVECTOR_HPP
#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <array>
#include <cstddef>

namespace Opm {

/**
 * \brief A simple vector-like class that is usable in GPU kernels.
 * 
 * This will behave as a vector, but with a fixed maximum size.
 * 
 * \note This is not suitable for large data sets, as it uses a fixed-size array internally.
 */
template <typename T, std::size_t max_size>
class GpuFriendlyVector {
public:
    OPM_HOST_DEVICE GpuFriendlyVector(std::initializer_list<T> list) {
        OPM_ERROR_IF(list.size() > max_size, "Initializer list size exceeds maximum size");
        for (auto& element : list) {
            push_back(element);
        }
    }

    OPM_HOST_DEVICE GpuFriendlyVector() {}

    OPM_HOST_DEVICE void push_back(const T& value) {
        if (size_ < max_size) {
            data_[size_] = value;
            ++size_;
        } else {
            OPM_THROW(std::runtime_error, "Exceeded maximum size");
        }
    }

    OPM_HOST_DEVICE T& operator[](std::size_t index) {
        return data_[index];
    }

    OPM_HOST_DEVICE const T& operator[](std::size_t index) const {
        return data_[index];
    }

    OPM_HOST_DEVICE size_t size() const {
        return size_;
    }
    OPM_HOST_DEVICE void clear() {
        size_ = 0;
    }
    OPM_HOST_DEVICE bool empty() const {
        return size_ == 0;
    }
    
    OPM_HOST_DEVICE void resize(std::size_t new_size) {
        if (new_size > max_size) {
            OPM_THROW(std::runtime_error, "New size exceeds maximum size");
        }
        size_ = new_size;
    }
    
private:
    std::array<T, max_size> data_;
    std::size_t size_ = 0u;
};
} // namespace Opm
#endif // OPM_COMMON_UTILITY_GPUFRIENDLYVECTOR_HPP
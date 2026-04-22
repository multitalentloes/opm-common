// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright TODO ADD YEAR AND NAME OF AUTHOR

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
/*!
 * \file
 *
 * GPU-portable counterpart of \c Opm::EclSpecrockLawParams.
 *
 * The CPU \c EclSpecrockLawParams holds a \c Tabulated1DFunction backed by
 * \c std::vector and is therefore not trivially copyable to the device.
 * This class stores the same temperature/internal-energy sample arrays in
 * a templated \c Storage container (defaulting to a CPU vector) so it can
 * be instantiated as a \c VectorWithDefaultAllocator-backed CPU object,
 * a \c GpuBuffer-backed owning GPU object and a \c GpuView-backed
 * non-owning GPU object usable from a kernel.
 */
#ifndef OPM_GPU_ECL_SPECROCK_LAW_PARAMS_HPP
#define OPM_GPU_ECL_SPECROCK_LAW_PARAMS_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/VectorWithDefaultAllocator.hpp>
#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/common/utility/gpuistl_if_available.hpp>

#include <opm/material/common/EnsureFinalized.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>

namespace Opm {

template <class ScalarT, template <class> class Storage = ::Opm::VectorWithDefaultAllocator>
class GpuEclSpecrockLawParams;

} // namespace Opm

#if HAVE_CUDA
namespace Opm::gpuistl {

template <class ScalarT>
::Opm::GpuEclSpecrockLawParams<ScalarT, GpuBuffer>
copy_to_gpu(const ::Opm::GpuEclSpecrockLawParams<ScalarT>& cpu);

template <class ScalarT, template <class> class ContainerT>
::Opm::GpuEclSpecrockLawParams<ScalarT, GpuView>
make_view(::Opm::GpuEclSpecrockLawParams<ScalarT, ContainerT>& gpuBuffers);

} // namespace Opm::gpuistl
#endif // HAVE_CUDA

namespace Opm {

/*!
 * \brief GPU-portable parameter object for the SPECROCK solid energy law.
 *
 * Stores the temperature-vs-volumetric-internal-energy table per cell or
 * per SATNUM region. The \c eval member performs piecewise-linear
 * interpolation matching \c Tabulated1DFunction::eval, with bisection over
 * the sample arrays.
 *
 * \tparam ScalarT  Floating point type used for the sample tables.
 * \tparam Storage  Template-template container holding the sample arrays.
 *                  Defaults to \c VectorWithDefaultAllocator (CPU vectors).
 */
template <class ScalarT, template <class> class Storage>
class GpuEclSpecrockLawParams : public EnsureFinalized
{
public:
    using Scalar = ScalarT;
    using ValueVector = Storage<Scalar>;

    OPM_HOST_DEVICE GpuEclSpecrockLawParams() = default;

    OPM_HOST_DEVICE GpuEclSpecrockLawParams(ValueVector temperatureSamples,
                                            ValueVector internalEnergySamples)
        : temperatureSamples_(std::move(temperatureSamples))
        , internalEnergySamples_(std::move(internalEnergySamples))
    {
        EnsureFinalized::finalize();
    }

    /*!
     * \brief Set the sample tables from any container with \c size() and
     *        \c operator[]. Marks the object as finalized.
     *
     * Available only on the CPU instantiation since GPU storage types are
     * not constructible from arbitrary host containers.
     */
    template <class ContainerT,
              class StorageT = Storage<Scalar>,
              std::enable_if_t<std::is_same_v<StorageT, ::Opm::VectorWithDefaultAllocator<Scalar>>,
                               int> = 0>
    void setSamples(const ContainerT& temperature, const ContainerT& internalEnergy)
    {
        if (temperature.size() != internalEnergy.size()) {
            OPM_THROW(std::invalid_argument,
                      "GpuEclSpecrockLawParams: temperature and internal-energy arrays must have "
                      "matching sizes");
        }
        const std::size_t n = temperature.size();
        temperatureSamples_.resize(n);
        internalEnergySamples_.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            temperatureSamples_[i] = static_cast<Scalar>(temperature[i]);
            internalEnergySamples_[i] = static_cast<Scalar>(internalEnergy[i]);
        }
        EnsureFinalized::finalize();
    }

    OPM_HOST_DEVICE std::size_t numSamples() const
    {
        return temperatureSamples_.size();
    }

    /*!
     * \brief Linearly interpolate the volumetric internal energy at a
     *        given temperature. The sample table is assumed sorted in
     *        ascending order; values outside the range are extrapolated
     *        linearly using the first/last segment (matching the CPU
     *        \c EclSpecrockLaw which calls \c eval(T, extrapolate=true)).
     */
    template <class Evaluation>
    OPM_HOST_DEVICE Evaluation eval(const Evaluation& x) const
    {
        const std::size_t n = temperatureSamples_.size();
        // n >= 2 by construction (SPECROCK tables always have >= 2 rows).
        std::size_t segIdx = 0;
        if (x <= temperatureSamples_[1]) {
            segIdx = 0;
        } else if (x >= temperatureSamples_[n - 2]) {
            segIdx = n - 2;
        } else {
            std::size_t lo = 1;
            std::size_t hi = n - 2;
            while (lo + 1 < hi) {
                const std::size_t mid = (lo + hi) / 2;
                if (x < temperatureSamples_[mid]) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            segIdx = lo;
        }
        const Scalar x0 = temperatureSamples_[segIdx];
        const Scalar x1 = temperatureSamples_[segIdx + 1];
        const Scalar y0 = internalEnergySamples_[segIdx];
        const Scalar y1 = internalEnergySamples_[segIdx + 1];
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
    }

    OPM_HOST_DEVICE const ValueVector& temperatureSamples() const
    {
        return temperatureSamples_;
    }

    OPM_HOST_DEVICE const ValueVector& internalEnergySamples() const
    {
        return internalEnergySamples_;
    }

    ValueVector& temperatureSamplesMutable()
    {
        return temperatureSamples_;
    }

    ValueVector& internalEnergySamplesMutable()
    {
        return internalEnergySamples_;
    }

private:
    ValueVector temperatureSamples_ {};
    ValueVector internalEnergySamples_ {};
};

} // namespace Opm

#if HAVE_CUDA
namespace Opm::gpuistl {

template <class ScalarT>
::Opm::GpuEclSpecrockLawParams<ScalarT, GpuBuffer>
copy_to_gpu(const ::Opm::GpuEclSpecrockLawParams<ScalarT>& cpu)
{
    return ::Opm::GpuEclSpecrockLawParams<ScalarT, GpuBuffer>(
        GpuBuffer<ScalarT>(cpu.temperatureSamples()),
        GpuBuffer<ScalarT>(cpu.internalEnergySamples()));
}

template <class ScalarT, template <class> class ContainerT>
::Opm::GpuEclSpecrockLawParams<ScalarT, GpuView>
make_view(::Opm::GpuEclSpecrockLawParams<ScalarT, ContainerT>& gpuBuffers)
{
    auto tView = make_view<ScalarT>(gpuBuffers.temperatureSamplesMutable());
    auto eView = make_view<ScalarT>(gpuBuffers.internalEnergySamplesMutable());
    return ::Opm::GpuEclSpecrockLawParams<ScalarT, GpuView>(tView, eView);
}

} // namespace Opm::gpuistl
#endif // HAVE_CUDA

#endif // OPM_GPU_ECL_SPECROCK_LAW_PARAMS_HPP

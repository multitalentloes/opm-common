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
 * GPU-portable counterpart of \c Opm::EclSpecrockLaw. Mirrors the static
 * \c solidInternalEnergy interface used by
 * \c BlackOilEnergyIntensiveQuantities::updateEnergyQuantities_, but is
 * decorated for execution from the device and operates on
 * \c Opm::GpuEclSpecrockLawParams.
 */
#ifndef OPM_GPU_ECL_SPECROCK_LAW_HPP
#define OPM_GPU_ECL_SPECROCK_LAW_HPP

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/thermal/GpuEclSpecrockLawParams.hpp>

namespace Opm {

/*!
 * \brief GPU-friendly volumetric internal-energy law for the rock,
 *        matching ECL SPECROCK semantics.
 */
template <class ScalarT, class ParamsT = GpuEclSpecrockLawParams<ScalarT>>
class GpuEclSpecrockLaw
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    OPM_HOST_DEVICE static Evaluation solidInternalEnergy(const Params& params,
                                                          const FluidState& fluidState)
    {
        const auto& temperature = fluidState.temperature(/*phaseIdx=*/0);
        return params.template eval<Evaluation>(temperature);
    }
};

} // namespace Opm

#endif // OPM_GPU_ECL_SPECROCK_LAW_HPP

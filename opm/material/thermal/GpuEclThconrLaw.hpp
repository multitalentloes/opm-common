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
 * GPU-portable counterpart of \c Opm::EclThconrLaw. The CPU
 * \c EclThconrLawParams is already trivially copyable (two scalars only)
 * and is reused here directly. The static \c thermalConductivity entry
 * point mirrors the CPU class but is \c OPM_HOST_DEVICE-decorated so it
 * can be invoked from a kernel.
 */
#ifndef OPM_GPU_ECL_THCONR_LAW_HPP
#define OPM_GPU_ECL_THCONR_LAW_HPP

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/densead/Math.hpp>
#include <opm/material/thermal/EclThconrLawParams.hpp>

namespace Opm {

/*!
 * \brief GPU-friendly total-thermal-conductivity law mirroring ECL THCONR
 *        + THCONSF.
 */
template <class ScalarT,
          class FluidSystem,
          class ParamsT = EclThconrLawParams<ScalarT>>
class GpuEclThconrLaw
{
public:
    using Params = ParamsT;
    using Scalar = typename Params::Scalar;

    template <class FluidState, class Evaluation = typename FluidState::ValueType>
    OPM_HOST_DEVICE static Evaluation thermalConductivity(const Params& params,
                                                          const FluidState& fluidState)
    {
        const Scalar lambdaRef = params.referenceTotalThermalConductivity();
        constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
        if (fluidState.fluidSystem().phaseIsActive(gasPhaseIdx)) {
            const Scalar alpha = params.dTotalThermalConductivity_dSg();
            const Evaluation& Sg = decay<Evaluation>(fluidState.saturation(gasPhaseIdx));
            return lambdaRef * (Scalar(1) - alpha * Sg);
        }
        return Evaluation(lambdaRef);
    }
};

} // namespace Opm

#endif // OPM_GPU_ECL_THCONR_LAW_HPP

// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
/*!
 * \file
 * \copydoc Opm::EclTwoPhaseMaterialParams
 */
#ifndef OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP
#define OPM_ECL_TWO_PHASE_MATERIAL_PARAMS_HPP

#include <memory>

#include <opm/material/common/EnsureFinalized.hpp>

#include <opm/common/utility/gpuDecorators.hpp>

namespace Opm {

enum class EclTwoPhaseApproach {
    GasOil,
    OilWater,
    GasWater
};

/*!
 * \brief Implementation for the parameters required by the material law for two-phase
 *        simulations.
 *
 * Essentially, this class just stores the two parameter objects for
 * the twophase capillary pressure laws.
 */
template<class Traits,
    class GasOilParamsT,
    class OilWaterParamsT,
    class GasWaterParamsT,
    template<typename> class SmartPointer = std::shared_ptr>
class EclTwoPhaseMaterialParams : public EnsureFinalized
{
    using Scalar = typename Traits::Scalar;
    enum { numPhases = 3 };
public:
    using EnsureFinalized :: finalize;

    using GasOilParams = GasOilParamsT;
    using OilWaterParams = OilWaterParamsT;
    using GasWaterParams = GasWaterParamsT;

    /*!
     * \brief The default constructor.
     */
    OPM_HOST_DEVICE EclTwoPhaseMaterialParams()
    {
    }

    OPM_HOST_DEVICE void setApproach(EclTwoPhaseApproach newApproach)
    { approach_ = newApproach; }

    OPM_HOST_DEVICE EclTwoPhaseApproach approach() const
    { return approach_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    OPM_HOST_DEVICE const GasOilParams& gasOilParams() const
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief The parameter object for the gas-oil twophase law.
     */
    OPM_HOST_DEVICE GasOilParams& gasOilParams()
    { EnsureFinalized::check(); return *gasOilParams_; }

    /*!
     * \brief Set the parameter object for the gas-oil twophase law.
     */
    OPM_HOST_DEVICE void setGasOilParams(SmartPointer<GasOilParams> val)
    { gasOilParams_ = val; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    OPM_HOST_DEVICE const OilWaterParams& oilWaterParams() const
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief The parameter object for the oil-water twophase law.
     */
    OPM_HOST_DEVICE OilWaterParams& oilWaterParams()
    { EnsureFinalized::check(); return *oilWaterParams_; }

    /*!
     * \brief Set the parameter object for the oil-water twophase law.
     */
    OPM_HOST_DEVICE void setOilWaterParams(SmartPointer<OilWaterParams> val)
    { oilWaterParams_ = val; }

  /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    OPM_HOST_DEVICE const GasWaterParams& gasWaterParams() const
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief The parameter object for the gas-water twophase law.
     */
    OPM_HOST_DEVICE GasWaterParams& gasWaterParams()
    { EnsureFinalized::check(); return *gasWaterParams_; }

    /*!
     * \brief Set the parameter object for the gas-water twophase law.
     */
    OPM_HOST_DEVICE void setGasWaterParams(SmartPointer<GasWaterParams> val)
    { gasWaterParams_ = val; }

    template<class Serializer>
    OPM_HOST_DEVICE void serializeOp(Serializer& serializer)
    {
        // This is for restart serialization.
        // Only dynamic state in the parameters need to be stored.
        serializer(*gasOilParams_);
        serializer(*oilWaterParams_);
        serializer(*gasWaterParams_);
    }

    OPM_HOST_DEVICE void setSwl(Scalar) {}

    OPM_HOST_DEVICE SmartPointer<GasOilParams>& gasOilParamsPtr() { return gasOilParams_; }
    OPM_HOST_DEVICE SmartPointer<OilWaterParams>& oilWaterParamsPtr() { return oilWaterParams_; }
    OPM_HOST_DEVICE SmartPointer<GasWaterParams>& gasWaterParamsPtr() { return gasWaterParams_; }

    OPM_HOST_DEVICE const SmartPointer<GasOilParams>& gasOilParamsPtr() const { return gasOilParams_; }
    OPM_HOST_DEVICE const SmartPointer<OilWaterParams>& oilWaterParamsPtr() const { return oilWaterParams_; }
    OPM_HOST_DEVICE const SmartPointer<GasWaterParams>& gasWaterParamsPtr() const { return gasWaterParams_; }

private:
    EclTwoPhaseApproach approach_{EclTwoPhaseApproach::GasOil};

    SmartPointer<GasOilParams> gasOilParams_;
    SmartPointer<OilWaterParams> oilWaterParams_;
    SmartPointer<GasWaterParams> gasWaterParams_;
};

namespace gpuistl {
    template<class Traits, class GPUContainerDouble, class GPUContainerScalar, class Scalar>
    EclTwoPhaseMaterialParams<Traits, GPUContainerDouble, GPUContainerScalar>
    copy_to_gpu(EclTwoPhaseMaterialParams<Traits, std::vector<double>, std::vector<Scalar>>& materialParams)
    {
        using GasOilParamsGPU = typename EclTwoPhaseMaterialParams<Traits, GPUContainerDouble, GPUContainerScalar>::GasOilParams;
        using OilWaterParamsGPU = typename EclTwoPhaseMaterialParams<Traits, GPUContainerDouble, GPUContainerScalar>::OilWaterParams;
        using GasWaterParamsGPU = typename EclTwoPhaseMaterialParams<Traits, GPUContainerDouble, GPUContainerScalar>::GasWaterParams;

        return EclTwoPhaseMaterialParams<Traits, GPUContainerDouble, GPUContainerScalar>(
            materialParams.approach(),
            copy_to_gpu<GPUContainerScalar, GasOilParamsGPU>(*materialParams.gasOilParams()),
            copy_to_gpu<GPUContainerScalar, OilWaterParamsGPU>(*materialParams.oilWaterParams()),
            copy_to_gpu<GPUContainerScalar, GasWaterParamsGPU>(*materialParams.gasWaterParams())
        );
    }

    template<template<class> class ViewPtr, class ViewDouble, class ViewScalar, class GPUContainerDouble, class GPUContainerScalar, class Traits>
    EclTwoPhaseMaterialParams<Traits, ViewDouble, ViewScalar, ViewPtr>
    make_view(EclTwoPhaseMaterialParams<Traits, GPUContainerDouble, GPUContainerScalar, std::unique_ptr>& materialParams)
    {
        using GasOilParamsView = typename EclTwoPhaseMaterialParams<Traits, ViewDouble, ViewScalar, ViewPtr>::GasOilParams;
        using OilWaterParamsView = typename EclTwoPhaseMaterialParams<Traits, ViewDouble, ViewScalar, ViewPtr>::OilWaterParams;
        using GasWaterParamsView = typename EclTwoPhaseMaterialParams<Traits, ViewDouble, ViewScalar, ViewPtr>::GasWaterParams;

        return EclTwoPhaseMaterialParams<Traits, ViewDouble, ViewScalar, ViewPtr>(
            materialParams.approach(),
            make_view<ViewScalar, GasOilParamsView>(*materialParams.gasOilParams()),
            make_view<ViewScalar, OilWaterParamsView>(*materialParams.oilWaterParams()),
            make_view<ViewScalar, GasWaterParamsView>(*materialParams.gasWaterParams())
        );
    }
}

} // namespace Opm

#endif

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
    template<class> class PtrType = std::shared_ptr>
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
    OPM_HOST_DEVICE void setGasOilParams(PtrType<GasOilParams> val)
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
    OPM_HOST_DEVICE void setOilWaterParams(PtrType<OilWaterParams> val)
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
    OPM_HOST_DEVICE void setGasWaterParams(PtrType<GasWaterParams> val)
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

    OPM_HOST_DEVICE PtrType<GasOilParams>& gasOilParamsPtr() { return gasOilParams_; }
    OPM_HOST_DEVICE PtrType<OilWaterParams>& oilWaterParamsPtr() { return oilWaterParams_; }
    OPM_HOST_DEVICE PtrType<GasWaterParams>& gasWaterParamsPtr() { return gasWaterParams_; }

    OPM_HOST_DEVICE const PtrType<GasOilParams>& gasOilParamsPtr() const { return gasOilParams_; }
    OPM_HOST_DEVICE const PtrType<OilWaterParams>& oilWaterParamsPtr() const { return oilWaterParams_; }
    OPM_HOST_DEVICE const PtrType<GasWaterParams>& gasWaterParamsPtr() const { return gasWaterParams_; }

private:
    EclTwoPhaseApproach approach_{EclTwoPhaseApproach::GasOil};

    PtrType<GasOilParams> gasOilParams_;
    PtrType<OilWaterParams> oilWaterParams_;
    PtrType<GasWaterParams> gasWaterParams_;
};

namespace gpuistl
{
    template<class ScalarGpuBuffer,
             class NewGasOilParamsT,
             class NewOilWaterParamsT,
             class NewGasWaterParamsT,
             class Traits,
             class OldGasOilParamsT,
             class OldOilWaterParamsT,
             class OldGasWaterParamsT>
    EclTwoPhaseMaterialParams<Traits, NewGasOilParamsT, NewOilWaterParamsT, NewGasWaterParamsT>
    copy_to_gpu(const EclTwoPhaseMaterialParams<Traits, OldGasOilParamsT, OldOilWaterParamsT, OldGasWaterParamsT>& params)
    {
        // Maybe I will run into some proble on a twophase case where some of these do not exist?
        // copy interpolation tables to the GPU - right now assumed to be the piecewiselinear....params
        auto gasOilParams = gpuistl::copy_to_gpu<ScalarGpuBuffer>(*params.gasOilParamsPtr());
        auto oilWaterParams = gpuistl::copy_to_gpu<ScalarGpuBuffer>(*params.oilWaterParamsPtr());
        auto gasWaterParams = gpuistl::copy_to_gpu<ScalarGpuBuffer>(*params.gasWaterParamsPtr());
        // Wrap the copied parameters in a shared_ptr
        auto gasOilParamsPtr = std::make_shared<NewGasOilParamsT>(gasOilParams);
        auto oilWaterParamsPtr = std::make_shared<NewOilWaterParamsT>(oilWaterParams);
        auto gasWaterParamsPtr = std::make_shared<NewGasWaterParamsT>(gasWaterParams);

        // Create the new EclTwoPhaseMaterialParams object
        auto gpuBufBasedEclTwoPhasedMaterialParams =
            EclTwoPhaseMaterialParams<Traits, NewGasOilParamsT, NewOilWaterParamsT, NewGasWaterParamsT>();

        gpuBufBasedEclTwoPhasedMaterialParams.setApproach(params.approach());
        gpuBufBasedEclTwoPhasedMaterialParams.setGasOilParams(gasOilParamsPtr);
        gpuBufBasedEclTwoPhasedMaterialParams.setOilWaterParams(oilWaterParamsPtr);
        gpuBufBasedEclTwoPhasedMaterialParams.setGasWaterParams(gasWaterParamsPtr);

        return gpuBufBasedEclTwoPhasedMaterialParams;
    }

    template<class ViewType,
             class NewGasOilParamsT,
            class NewOilWaterParamsT,
            class NewGasWaterParamsT,
            template<class> class PtrType,
            class Traits,
            class OldGasOilParamsT,
            class OldOilWaterParamsT,
            class OldGasWaterParamsT>
    EclTwoPhaseMaterialParams<Traits, NewGasOilParamsT, NewOilWaterParamsT, NewGasWaterParamsT, PtrType>
    make_view(const EclTwoPhaseMaterialParams<Traits, OldGasOilParamsT, OldOilWaterParamsT, OldGasWaterParamsT>& params)
    {
        // Maybe I will run into some proble on a twophase case where some of these do not exist?
        // copy interpolation tables to the GPU - right now assumed to be the piecewiselinear....params
        auto gasOilParams = gpuistl::make_view<ViewType>(*params.gasOilParamsPtr());
        auto oilWaterParams = gpuistl::make_view<ViewType>(*params.oilWaterParamsPtr());
        auto gasWaterParams = gpuistl::make_view<ViewType>(*params.gasWaterParamsPtr());

        // Not immediately convinced this is correct
        auto gasOilParamsPtr = PtrType<NewGasOilParamsT>(gasOilParams);
        auto oilWaterParamsPtr = PtrType<NewOilWaterParamsT>(oilWaterParams);
        auto gasWaterParamsPtr = PtrType<NewGasWaterParamsT>(gasWaterParams);

        // Create the new EclTwoPhaseMaterialParams object
        auto gpuBufBasedEclTwoPhasedMaterialParams =
            EclTwoPhaseMaterialParams<Traits, NewGasOilParamsT, NewOilWaterParamsT, NewGasWaterParamsT, PtrType>();

        gpuBufBasedEclTwoPhasedMaterialParams.setApproach(params.approach());
        gpuBufBasedEclTwoPhasedMaterialParams.setGasOilParams(gasOilParamsPtr);
        gpuBufBasedEclTwoPhasedMaterialParams.setOilWaterParams(oilWaterParamsPtr);
        gpuBufBasedEclTwoPhasedMaterialParams.setGasWaterParams(gasWaterParamsPtr);

        return gpuBufBasedEclTwoPhasedMaterialParams;
    }
}

} // namespace Opm

#endif

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
 *
 * \copydoc Opm::BlackOilFluidState
 */
#ifndef OPM_BLACK_OIL_FLUID_STATE_HH
#define OPM_BLACK_OIL_FLUID_STATE_HH

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/common/HasMemberGeneratorMacros.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/ConditionalStorage.hpp>

#include <opm/common/utility/gpuDecorators.hpp>

namespace Opm {

// TODO: do I have to do anything for GPU to accomodate this macro?
OPM_GENERATE_HAS_MEMBER(pvtRegionIndex, ) // Creates 'HasMember_pvtRegionIndex<T>'.

template <class FluidState>
OPM_HOST_DEVICE unsigned getPvtRegionIndex_(typename std::enable_if<HasMember_pvtRegionIndex<FluidState>::value,
                                                    const FluidState&>::type fluidState)
{ return fluidState.pvtRegionIndex(); }

template <class FluidState>
OPM_HOST_DEVICE unsigned getPvtRegionIndex_(typename std::enable_if<!HasMember_pvtRegionIndex<FluidState>::value,
                                                    const FluidState&>::type)
{ return 0; }

OPM_GENERATE_HAS_MEMBER(invB, /*phaseIdx=*/0) // Creates 'HasMember_invB<T>'.

template <class FluidSystem, class FluidState, class LhsEval>
OPM_HOST_DEVICE auto getInvB_(typename std::enable_if<HasMember_invB<FluidState>::value,
                                      const FluidState&>::type fluidState,
              unsigned phaseIdx,
              unsigned)
    -> decltype(decay<LhsEval>(fluidState.invB(phaseIdx)))
{ return decay<LhsEval>(fluidState.invB(phaseIdx)); }

template <class FluidSystem, class FluidState, class LhsEval>
OPM_HOST_DEVICE LhsEval getInvB_(typename std::enable_if<!HasMember_invB<FluidState>::value,
                                         const FluidState&>::type fluidState,
                 unsigned phaseIdx,
                 unsigned pvtRegionIdx)
{
    const auto& rho = fluidState.density(phaseIdx);
    const auto& Xsolvent =
        fluidState.massFraction(phaseIdx, FluidSystem::solventComponentIndex(phaseIdx));

    return
        decay<LhsEval>(rho)
        *decay<LhsEval>(Xsolvent)
        /FluidSystem::referenceDensity(phaseIdx, pvtRegionIdx);
}

OPM_GENERATE_HAS_MEMBER(saltConcentration, ) // Creates 'HasMember_saltConcentration<T>'.

template <class FluidState>
OPM_HOST_DEVICE auto getSaltConcentration_(typename std::enable_if<HasMember_saltConcentration<FluidState>::value,
                                                    const FluidState&>::type fluidState)
{ return fluidState.saltConcentration(); }

template <class FluidState>
OPM_HOST_DEVICE auto getSaltConcentration_(typename std::enable_if<!HasMember_saltConcentration<FluidState>::value,
                                                    const FluidState&>::type)
{ return 0.0; }

OPM_GENERATE_HAS_MEMBER(saltSaturation, ) // Creates 'HasMember_saltSaturation<T>'.

template <class FluidState>
OPM_HOST_DEVICE auto getSaltSaturation_(typename std::enable_if<HasMember_saltSaturation<FluidState>::value,
                                                    const FluidState&>::type fluidState)
{ return fluidState.saltSaturation(); }


template <class FluidState>
OPM_HOST_DEVICE auto getSaltSaturation_(typename std::enable_if<!HasMember_saltSaturation<FluidState>::value,
                                                    const FluidState&>::type)
{ return 0.0; }

/*!
 * \brief Implements a "tailor-made" fluid state class for the black-oil model.
 *
 * I.e., it uses exactly the same quantities which are used by the ECL blackoil
 * model. Further quantities are computed "on the fly" and are accessing them is thus
 * relatively slow.
 */
template <class ScalarT,
          class FluidSystem,
          bool enableTemperature = false,
          bool enableEnergy = false,
          bool enableDissolution = true,
          bool enableVapwat = false,
          bool enableBrine = false,
          bool enableSaltPrecipitation = false,
          bool enableDissolutionInWater = false,
          unsigned numStoragePhases = FluidSystem::numPhases>
class BlackOilFluidState
{
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };

    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

public:
    using Scalar = ScalarT;
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };

    /*!
     * \brief Make sure that all attributes are defined.
     *
     * This method does not do anything if the program is not run
     * under valgrind. If it is, then valgrind will print an error
     * message if some attributes of the object have not been properly
     * defined.
     */
    OPM_HOST_DEVICE void checkDefined() const
    {
#ifndef NDEBUG
#if !OPM_IS_INSIDE_DEVICE_FUNCTION
        Valgrind::CheckDefined(pvtRegionIdx_);

        for (unsigned storagePhaseIdx = 0; storagePhaseIdx < numStoragePhases; ++ storagePhaseIdx) {
            Valgrind::CheckDefined(saturation_[storagePhaseIdx]);
            Valgrind::CheckDefined(pressure_[storagePhaseIdx]);
            Valgrind::CheckDefined(density_[storagePhaseIdx]);
            Valgrind::CheckDefined(invB_[storagePhaseIdx]);

            if constexpr (enableEnergy)
                Valgrind::CheckDefined((*enthalpy_)[storagePhaseIdx]);
        }

        if constexpr (enableDissolution) {
            Valgrind::CheckDefined(*Rs_);
            Valgrind::CheckDefined(*Rv_);
        }

        if constexpr (enableVapwat) {
            Valgrind::CheckDefined(*Rvw_);
        }

        if constexpr (enableDissolutionInWater) {
            Valgrind::CheckDefined(*Rsw_);
        }

        if constexpr (enableBrine) {
            Valgrind::CheckDefined(*saltConcentration_);
        }

        if constexpr (enableSaltPrecipitation) {
            Valgrind::CheckDefined(*saltSaturation_);
        }

        if constexpr (enableTemperature || enableEnergy)
            Valgrind::CheckDefined(*temperature_);
#endif // !OPM_IS_INSIDE_DEVICE_FUNCTION
#endif // NDEBUG
    }

    /*!
     * \brief Retrieve all parameters from an arbitrary fluid
     *        state.
     */
    template <class FluidState>
    OPM_HOST_DEVICE void assign(const FluidState& fs)
    {
        if constexpr (enableTemperature || enableEnergy)
            setTemperature(fs.temperature(/*phaseIdx=*/0));

        unsigned pvtRegionIdx = getPvtRegionIndex_<FluidState>(fs);
        setPvtRegionIndex(pvtRegionIdx);

        if constexpr (enableDissolution) {
            setRs(BlackOil::getRs_<FluidSystem, FluidState, Scalar>(fs, pvtRegionIdx));
            setRv(BlackOil::getRv_<FluidSystem, FluidState, Scalar>(fs, pvtRegionIdx));
        }
        if constexpr (enableVapwat) {
            setRvw(BlackOil::getRvw_<FluidSystem, FluidState, Scalar>(fs, pvtRegionIdx));
        }
        if constexpr (enableDissolutionInWater) {
            setRsw(BlackOil::getRsw_<FluidSystem, FluidState, Scalar>(fs, pvtRegionIdx));
        }
        if constexpr (enableBrine){
            setSaltConcentration(BlackOil::getSaltConcentration_<FluidSystem, FluidState, Scalar>(fs, pvtRegionIdx));
        }
        if constexpr (enableSaltPrecipitation){
            setSaltSaturation(BlackOil::getSaltSaturation_<FluidSystem, FluidState, Scalar>(fs, pvtRegionIdx));
        }
        for (unsigned storagePhaseIdx = 0; storagePhaseIdx < numStoragePhases; ++storagePhaseIdx) {
            unsigned phaseIdx = storageToCanonicalPhaseIndex_(storagePhaseIdx);
            setSaturation(phaseIdx, fs.saturation(phaseIdx));
            setPressure(phaseIdx, fs.pressure(phaseIdx));
            setDensity(phaseIdx, fs.density(phaseIdx));

            if constexpr (enableEnergy)
                setEnthalpy(phaseIdx, fs.enthalpy(phaseIdx));

            setInvB(phaseIdx, getInvB_<FluidSystem, FluidState, Scalar>(fs, phaseIdx, pvtRegionIdx));
        }
    }

    /*!
     * \brief Set the index of the fluid region
     *
     * This determines which tables are used to compute the quantities that are computed
     * on the fly.
     */
    OPM_HOST_DEVICE void setPvtRegionIndex(unsigned newPvtRegionIdx)
    { pvtRegionIdx_ = static_cast<unsigned short>(newPvtRegionIdx); }

    /*!
     * \brief Set the pressure of a fluid phase [-].
     */
    OPM_HOST_DEVICE void setPressure(unsigned phaseIdx, const Scalar& p)
    { pressure_[canonicalToStoragePhaseIndex_(phaseIdx)] = p; }

    /*!
     * \brief Set the saturation of a fluid phase [-].
     */
    OPM_HOST_DEVICE void setSaturation(unsigned phaseIdx, const Scalar& S)
    { saturation_[canonicalToStoragePhaseIndex_(phaseIdx)] = S; }

    /*!
     * \brief Set the capillary pressure of a fluid phase [-].
     */
    OPM_HOST_DEVICE void setPc(unsigned phaseIdx, const Scalar& pc)
    { pc_[canonicalToStoragePhaseIndex_(phaseIdx)] = pc; }

    /*!
     * \brief Set the total saturation used for sequential methods
     */
    OPM_HOST_DEVICE void setTotalSaturation(const Scalar& value)
    {
        totalSaturation_ = value;
    }

    /*!
     * \brief Set the temperature [K]
     *
     * If neither the enableTemperature nor the enableEnergy template arguments are set
     * to true, this method will throw an exception!
     */
    OPM_HOST_DEVICE void setTemperature(const Scalar& value)
    {
        assert(enableTemperature || enableEnergy);

        (*temperature_) = value;
    }

    /*!
     * \brief Set the specific enthalpy [J/kg] of a given fluid phase.
     *
     * If the enableEnergy template argument is not set to true, this method will throw
     * an exception!
     */
    OPM_HOST_DEVICE void setEnthalpy(unsigned phaseIdx, const Scalar& value)
    {
        assert(enableTemperature || enableEnergy);

        (*enthalpy_)[canonicalToStoragePhaseIndex_(phaseIdx)] = value;
    }

    /*!
     * \ brief Set the inverse formation volume factor of a fluid phase
     */
    OPM_HOST_DEVICE void setInvB(unsigned phaseIdx, const Scalar& b)
    { invB_[canonicalToStoragePhaseIndex_(phaseIdx)] = b; }

    /*!
     * \ brief Set the density of a fluid phase
     */
    OPM_HOST_DEVICE void setDensity(unsigned phaseIdx, const Scalar& rho)
    { density_[canonicalToStoragePhaseIndex_(phaseIdx)] = rho; }

    /*!
     * \brief Set the gas dissolution factor [m^3/m^3] of the oil phase.
     *
     * This quantity is very specific to the black-oil model.
     */
    OPM_HOST_DEVICE void setRs(const Scalar& newRs)
    { *Rs_ = newRs; }

    /*!
     * \brief Set the oil vaporization factor [m^3/m^3] of the gas phase.
     *
     * This quantity is very specific to the black-oil model.
     */
    OPM_HOST_DEVICE void setRv(const Scalar& newRv)
    { *Rv_ = newRv; }

    /*!
     * \brief Set the water vaporization factor [m^3/m^3] of the gas phase.
     *
     * This quantity is very specific to the black-oil model.
     */
    OPM_HOST_DEVICE void setRvw(const Scalar& newRvw)
    { *Rvw_ = newRvw; }

    /*!
     * \brief Set the gas dissolution factor [m^3/m^3] of the water phase..
     *
     * This quantity is very specific to the black-oil model.
     */
    OPM_HOST_DEVICE void setRsw(const Scalar& newRsw)
    { *Rsw_ = newRsw; }

    /*!
     * \brief Set the salt concentration.
     */
    OPM_HOST_DEVICE void setSaltConcentration(const Scalar& newSaltConcentration)
    { *saltConcentration_ = newSaltConcentration; }

    /*!
     * \brief Set the solid salt saturation.
     */
    OPM_HOST_DEVICE void setSaltSaturation(const Scalar& newSaltSaturation)
    { *saltSaturation_ = newSaltSaturation; }

    /*!
     * \brief Return the pressure of a fluid phase [Pa]
     */
    OPM_HOST_DEVICE const Scalar& pressure(unsigned phaseIdx) const
    { return pressure_[canonicalToStoragePhaseIndex_(phaseIdx)]; }

    /*!
     * \brief Return the saturation of a fluid phase [-]
     */
    OPM_HOST_DEVICE const Scalar& saturation(unsigned phaseIdx) const
    { return saturation_[canonicalToStoragePhaseIndex_(phaseIdx)]; }

    /*!
     * \brief Return the capillary pressure of a fluid phase [-]
     */
    OPM_HOST_DEVICE const Scalar& pc(unsigned phaseIdx) const
    { return pc_[canonicalToStoragePhaseIndex_(phaseIdx)]; }

    /*!
     * \brief Return the total saturation needed for sequential
     */
    OPM_HOST_DEVICE const Scalar& totalSaturation() const
    {
        return totalSaturation_;
    }

    /*!
     * \brief Return the temperature [K]
     */
    OPM_HOST_DEVICE const Scalar& temperature(unsigned) const
    {
        if constexpr (enableTemperature || enableEnergy) {
            return *temperature_;
        } else {
            static Scalar tmp(FluidSystem::reservoirTemperature(pvtRegionIdx_));
            return tmp;
        }
    }

    /*!
     * \brief Return the inverse formation volume factor of a fluid phase [-].
     *
     * This factor expresses the change of density of a pure phase due to increased
     * pressure and temperature at reservoir conditions compared to surface conditions.
     */
    OPM_HOST_DEVICE const Scalar& invB(unsigned phaseIdx) const
    { return invB_[canonicalToStoragePhaseIndex_(phaseIdx)]; }

    /*!
     * \brief Return the gas dissolution factor of oil [m^3/m^3].
     *
     * I.e., the amount of gas which is present in the oil phase in terms of cubic meters
     * of gas at surface conditions per cubic meter of liquid oil at surface
     * conditions. This method is specific to the black-oil model.
     */
    OPM_HOST_DEVICE const Scalar& Rs() const
    {
        if constexpr (enableDissolution) {
            return *Rs_;
        } else {
            static Scalar null = 0.0;
            return null;
        }
    }

    /*!
     * \brief Return the oil vaporization factor of gas [m^3/m^3].
     *
     * I.e., the amount of oil which is present in the gas phase in terms of cubic meters
     * of liquid oil at surface conditions per cubic meter of gas at surface
     * conditions. This method is specific to the black-oil model.
     */
    OPM_HOST_DEVICE const Scalar& Rv() const
    {
        if constexpr (!enableDissolution) {
            static Scalar null = 0.0;
            return null;
        } else {
            return *Rv_;
        }
    }

    /*!
     * \brief Return the water vaporization factor of gas [m^3/m^3].
     *
     * I.e., the amount of water which is present in the gas phase in terms of cubic meters
     * of liquid water at surface conditions per cubic meter of gas at surface
     * conditions. This method is specific to the black-oil model.
     */
    OPM_HOST_DEVICE const Scalar& Rvw() const
    {
        if constexpr (enableVapwat) {
            return *Rvw_;
        } else {
            static Scalar null = 0.0;
            return null;
        }
    }

    /*!
     * \brief Return the gas dissolution factor of water [m^3/m^3].
     *
     * I.e., the amount of gas which is present in the water phase in terms of cubic meters
     * of gas at surface conditions per cubic meter of water at surface
     * conditions. This method is specific to the black-oil model.
     */
    OPM_HOST_DEVICE const Scalar& Rsw() const
    {
        if constexpr (enableDissolutionInWater) {
            return *Rsw_;
        } else {
            static Scalar null = 0.0;
            return null;
        }
    }

    /*!
     * \brief Return the concentration of salt in water
     */
    OPM_HOST_DEVICE const Scalar& saltConcentration() const
    {
        if constexpr (enableBrine) {
            return *saltConcentration_;
        } else {
            static Scalar null = 0.0;
            return null;
        }
    }

    /*!
     * \brief Return the saturation of solid salt
     */
    OPM_HOST_DEVICE const Scalar& saltSaturation() const
    {
        if constexpr (enableSaltPrecipitation) {
            return *saltSaturation_;
        } else {
            static Scalar null = 0.0;
            return null;
        }
    }

    /*!
     * \brief Return the PVT region where the current fluid state is assumed to be part of.
     *
     */
    OPM_HOST_DEVICE unsigned short pvtRegionIndex() const
    { return pvtRegionIdx_; }

    /*!
     * \brief Return the density [kg/m^3] of a given fluid phase.
      */
    OPM_HOST_DEVICE Scalar density(unsigned phaseIdx) const
    { return density_[canonicalToStoragePhaseIndex_(phaseIdx)]; }

    /*!
     * \brief Return the specific enthalpy [J/kg] of a given fluid phase.
     *
     * If the EnableEnergy property is not set to true, this method will throw an
     * exception!
     */
    OPM_HOST_DEVICE const Scalar& enthalpy(unsigned phaseIdx) const
    { return (*enthalpy_)[canonicalToStoragePhaseIndex_(phaseIdx)]; }

    /*!
     * \brief Return the specific internal energy [J/kg] of a given fluid phase.
     *
     * If the EnableEnergy property is not set to true, this method will throw an
     * exception!
     */
    OPM_HOST_DEVICE Scalar internalEnergy(unsigned phaseIdx) const
    {   auto energy = (*enthalpy_)[canonicalToStoragePhaseIndex_(phaseIdx)];
        if(!FluidSystem::enthalpyEqualEnergy()){
            energy -= pressure(phaseIdx)/density(phaseIdx);
        }
        return energy;
    }

    //////
    // slow methods
    //////

    /*!
     * \brief Return the molar density of a fluid phase [mol/m^3].
     */
    OPM_HOST_DEVICE Scalar molarDensity(unsigned phaseIdx) const
    {
        const auto& rho = density(phaseIdx);

        if (phaseIdx == waterPhaseIdx)
            return rho/FluidSystem::molarMass(waterCompIdx, pvtRegionIdx_);

        return
            rho*(moleFraction(phaseIdx, gasCompIdx)/FluidSystem::molarMass(gasCompIdx, pvtRegionIdx_)
                 + moleFraction(phaseIdx, oilCompIdx)/FluidSystem::molarMass(oilCompIdx, pvtRegionIdx_));

    }

    /*!
     * \brief Return the molar volume of a fluid phase [m^3/mol].
     *
     * This is equivalent to the inverse of the molar density.
     */
    OPM_HOST_DEVICE Scalar molarVolume(unsigned phaseIdx) const
    { return 1.0/molarDensity(phaseIdx); }

    /*!
     * \brief Return the dynamic viscosity of a fluid phase [Pa s].
     */
    OPM_HOST_DEVICE Scalar viscosity(unsigned phaseIdx) const
    { return FluidSystem::viscosity(*this, phaseIdx, pvtRegionIdx_); }

    /*!
     * \brief Return the mass fraction of a component in a fluid phase [-].
     */
    OPM_HOST_DEVICE Scalar massFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        switch (phaseIdx) {
        case waterPhaseIdx:
            if (compIdx == waterCompIdx)
                return 1.0;
            return 0.0;

        case oilPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return 1.0 - FluidSystem::convertRsToXoG(Rs(), pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return FluidSystem::convertRsToXoG(Rs(), pvtRegionIdx_);
            }
            break;

        case gasPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return FluidSystem::convertRvToXgO(Rv(), pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return 1.0 - FluidSystem::convertRvToXgO(Rv(), pvtRegionIdx_);
            }
            break;
        }

        throw std::logic_error("Invalid phase or component index!");
    }

    /*!
     * \brief Return the mole fraction of a component in a fluid phase [-].
     */
    OPM_HOST_DEVICE Scalar moleFraction(unsigned phaseIdx, unsigned compIdx) const
    {
        switch (phaseIdx) {
        case waterPhaseIdx:
            if (compIdx == waterCompIdx)
                return 1.0;
            return 0.0;

        case oilPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return 1.0 - FluidSystem::convertXoGToxoG(FluidSystem::convertRsToXoG(Rs(), pvtRegionIdx_),
                                                          pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return FluidSystem::convertXoGToxoG(FluidSystem::convertRsToXoG(Rs(), pvtRegionIdx_),
                                                    pvtRegionIdx_);
            }
            break;

        case gasPhaseIdx:
            if (compIdx == waterCompIdx)
                return 0.0;
            else if (compIdx == oilCompIdx)
                return FluidSystem::convertXgOToxgO(FluidSystem::convertRvToXgO(Rv(), pvtRegionIdx_),
                                                    pvtRegionIdx_);
            else {
                assert(compIdx == gasCompIdx);
                return 1.0 - FluidSystem::convertXgOToxgO(FluidSystem::convertRvToXgO(Rv(), pvtRegionIdx_),
                                                          pvtRegionIdx_);
            }
            break;
        }

        throw std::logic_error("Invalid phase or component index!");
    }

    /*!
     * \brief Return the partial molar density of a component in a fluid phase [mol / m^3].
     */
    OPM_HOST_DEVICE Scalar molarity(unsigned phaseIdx, unsigned compIdx) const
    { return moleFraction(phaseIdx, compIdx)*molarDensity(phaseIdx); }

    /*!
     * \brief Return the partial molar density of a fluid phase [kg / mol].
     */
    OPM_HOST_DEVICE Scalar averageMolarMass(unsigned phaseIdx) const
    {
        Scalar result(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++ compIdx)
            result += FluidSystem::molarMass(compIdx, pvtRegionIdx_)*moleFraction(phaseIdx, compIdx);
        return result;
    }

    /*!
     * \brief Return the fugacity coefficient of a component in a fluid phase [-].
     */
    OPM_HOST_DEVICE Scalar fugacityCoefficient(unsigned phaseIdx, unsigned compIdx) const
    { return FluidSystem::fugacityCoefficient(*this, phaseIdx, compIdx, pvtRegionIdx_); }

    /*!
     * \brief Return the fugacity of a component in a fluid phase [Pa].
     */
    OPM_HOST_DEVICE Scalar fugacity(unsigned phaseIdx, unsigned compIdx) const
    {
        return
            fugacityCoefficient(phaseIdx, compIdx)
            *moleFraction(phaseIdx, compIdx)
            *pressure(phaseIdx);
    }

private:
    OPM_HOST_DEVICE static unsigned storageToCanonicalPhaseIndex_(unsigned storagePhaseIdx)
    {
        if constexpr (numStoragePhases == 3)
            return storagePhaseIdx;
        else
            return FluidSystem::activeToCanonicalPhaseIdx(storagePhaseIdx);
    }

    OPM_HOST_DEVICE static unsigned canonicalToStoragePhaseIndex_(unsigned canonicalPhaseIdx)
    {
        if constexpr (numStoragePhases == 3)
            return canonicalPhaseIdx;
        else
            return FluidSystem::canonicalToActivePhaseIdx(canonicalPhaseIdx);
    }

    // All of this data is actually POD and we probably do not need copy_to_gpu and make_view

    ConditionalStorage<enableTemperature || enableEnergy, Scalar> temperature_{};
    ConditionalStorage<enableEnergy, std::array<Scalar, numStoragePhases> > enthalpy_{};
    Scalar totalSaturation_{};
    std::array<Scalar, numStoragePhases> pressure_{};
    std::array<Scalar, numStoragePhases> pc_{};
    std::array<Scalar, numStoragePhases> saturation_{};
    std::array<Scalar, numStoragePhases> invB_{};
    std::array<Scalar, numStoragePhases> density_{};
    ConditionalStorage<enableDissolution,Scalar> Rs_{};
    ConditionalStorage<enableDissolution, Scalar> Rv_{};
    ConditionalStorage<enableVapwat,Scalar> Rvw_{};
    ConditionalStorage<enableDissolutionInWater,Scalar> Rsw_{};
    ConditionalStorage<enableBrine, Scalar> saltConcentration_{};
    ConditionalStorage<enableSaltPrecipitation, Scalar> saltSaturation_{};

    unsigned short pvtRegionIdx_{};
};

} // namespace Opm

// namespace Opm::gpuistl {
//     template<class Scalar, class Params, class GPUContainer>
//     Co2GasPvt<Scalar, Params, GPUContainer>
//     copy_to_gpu(const Co2GasPvt<Scalar>& cpuCo2)
//     {
//         return Co2GasPvt<Scalar, Params, GPUContainer>(
//             copy_to_gpu<Scalar, std::vector<Scalar>, GPUContainer>(cpuCo2.getParams()),
//             GPUContainer(cpuCo2.getBrineReferenceDensity()),
//             GPUContainer(cpuCo2.getGasReferenceDensity()),
//             GPUContainer(cpuCo2.getSalinity()),
//             cpuCo2.getEnableEzrokhiDensity(),
//             cpuCo2.getEnableVaporization(),
//             cpuCo2.getActivityModel(),
//             cpuCo2.getGasType());
//     }

//     template <class ViewType, class OutputParams, class InputParams, class ContainerType, class Scalar>
//     Co2GasPvt<Scalar, OutputParams, ViewType>
//     make_view(const Co2GasPvt<Scalar, InputParams, ContainerType>& co2GasPvt)
//     {
//         using containedType = typename ContainerType::value_type;
//         using viewedTypeNoConst = typename std::remove_const_t<typename ViewType::value_type>;

//         static_assert(std::is_same_v<containedType, viewedTypeNoConst>);

//         ViewType newBrineReferenceDensity = make_view<viewedTypeNoConst>(co2GasPvt.getBrineReferenceDensity());
//         ViewType newGasReferenceDensity = make_view<viewedTypeNoConst>(co2GasPvt.getGasReferenceDensity());
//         ViewType newSalinity = make_view<viewedTypeNoConst>(co2GasPvt.getSalinity());

//         return Co2GasPvt<Scalar, OutputParams, ViewType>(
//             make_view<ViewType>(co2GasPvt.getParams()),
//             newBrineReferenceDensity,
//             newGasReferenceDensity,
//             newSalinity,
//             co2GasPvt.getEnableEzrokhiDensity(),
//             co2GasPvt.getEnableVaporization(),
//             co2GasPvt.getActivityModel(),
//             co2GasPvt.getGasType());
//     }
// }

#endif

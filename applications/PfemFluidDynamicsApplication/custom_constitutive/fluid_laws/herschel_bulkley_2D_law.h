//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Alessandro Franci, Ignasi de Pouplana
//  Collaborators:  Timur Tomas
//
//-------------------------------------------------------------
//

#if !defined(KRATOS_HERSCHEL_BULKLEY_LAW_2D_H_INCLUDED)
#define KRATOS_HERSCHEL_BULKLEY_LAW_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "fluid_constitutive_law.h"

namespace Kratos {
/**
 * Defines a 2D Herschel Bulkley non-Newtonian constitutive law
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 * 2) YIELD_SHEAR
 * 3) ADAPTIVE_EXPONENT
 * 4) FLOW_INDEX
 */
class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) HerschelBulkley2DLaw : public PfemFluidConstitutiveLaw {
   public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;

    /**
     * Counted pointer of HerschelBulkley2DLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(HerschelBulkley2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HerschelBulkley2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HerschelBulkley2DLaw(const HerschelBulkley2DLaw& rOther);

    /**
     * Destructor.
     */
    ~HerschelBulkley2DLaw() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * @return Working space dimension constitutive law
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @return Size of the strain vector (in Voigt notation) for the constitutive law
     */
    SizeType GetStrainSize();

    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo);

    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     */
    std::string Info() const override;

   protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /// Get the effective viscosity (in dynamic units -- Pa s) for the fluid.
    double GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const override;

    /// Get the effective density for the fluid.
    double GetEffectiveDensity(ConstitutiveLaw::Parameters& rParameters) const override;

    /// Get the effective yield shear for the fluid.
    double GetEffectiveYieldShear(ConstitutiveLaw::Parameters& rParameters) const;

    /// Get the effective dynamic viscosity for the fluid.
    double GetEffectiveDynamicViscosity(ConstitutiveLaw::Parameters& rParameters) const;

    /// Get the flow index for the fluid.
    double GetFlowIndex(ConstitutiveLaw::Parameters& rParameters) const;

    ///@}

   private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
    ///@}

};  // Class HerschelBulkley2DLaw

}  // namespace Kratos.

#endif  // KRATOS_HERSCHEL_BULKLEY_LAW_2D_H_INCLUDED  defined

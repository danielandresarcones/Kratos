//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PRECONDITIONER_FACTORY_H_INCLUDED )
#define  KRATOS_PRECONDITIONER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "includes/shared_pointers.h"
#include "linear_solvers/preconditioner.h"
#include "spaces/ublas_space.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class PreconditionerFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of preconditioners
 * @details Defines the base preconditioner factory
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 */
template< typename TSparseSpace, typename TLocalSpace>
class PreconditionerFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the custom class
    typedef PreconditionerFactory<TSparseSpace,TLocalSpace> FactoryType;

    /// The definition of the preconditioner
    typedef Preconditioner<TSparseSpace,TLocalSpace> PreconditionerType;

    /// Pointer definition of PreconditionerFactory
    KRATOS_CLASS_POINTER_DEFINITION(PreconditionerFactory );

    ///@}
    ///@name Life Cycle
    ///@{

    virtual ~PreconditionerFactory() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks if the linear solver is registered
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string& rSolverType) const
    {
        return KratosComponents< FactoryType >::Has( rSolverType );
    }

    /**
     * @brief This method creates a new solver
     * @return The pointer to the solver of interest
     */
    virtual typename PreconditionerType::Pointer Create(const std::string& rPreconditionerType)  const
    {
        // remove name of the application (if passed)
        // e.g. "LinearSolversApplication.sparse_lu" => "sparse_lu"
        const std::string raw_precond_name = rPreconditionerType.substr(rPreconditionerType.find(".") + 1);

        KRATOS_ERROR_IF_NOT(Has(raw_precond_name))
            << "Trying to construct a preconditioner with preconditioner_type:\n\""
            << raw_precond_name << "\" which does not exist.\n"
            << "The list of available options (for currently loaded applications) is:\n"
            << KratosComponents< FactoryType >() << std::endl;

        const auto& aux = KratosComponents< FactoryType >::Get(raw_precond_name);
        return aux.CreatePreconditioner();
    }
    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new solver
     * @return The pointer to the solver of interest
     */
    virtual typename PreconditionerType::Pointer CreatePreconditioner()  const
    {
        KRATOS_ERROR << "calling the base class PreconditionerFactory" << std::endl;
    }
};

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template< typename TSparseSpace, typename TLocalSpace>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PreconditionerFactory<TSparseSpace, TLocalSpace>& rThis)
{
    rOStream << "PreconditionerFactory" << std::endl;

    return rOStream;
}

///@}

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSparseSpaceType;

typedef PreconditionerFactory<SparseSpaceType, LocalSparseSpaceType> PreconditionerFactoryType;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<PreconditionerFactoryType>;

#ifdef KRATOS_REGISTER_PRECONDITIONER
#undef KRATOS_REGISTER_PRECONDITIONER
#endif
#define KRATOS_REGISTER_PRECONDITIONER(name, reference) \
    KratosComponents<PreconditionerFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_PRECONDITIONER_FACTORY_H_INCLUDED  defined

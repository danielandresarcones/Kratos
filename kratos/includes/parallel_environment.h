//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#ifndef KRATOS_PARALLEL_ENVIRONMENT_H_INCLUDED
#define KRATOS_PARALLEL_ENVIRONMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"
#include "includes/fill_communicator.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) EnvironmentManager
{
  public:
    typedef std::unique_ptr<EnvironmentManager> Pointer;

    EnvironmentManager(EnvironmentManager& rOther) = delete;

    virtual ~EnvironmentManager() = default;

    virtual bool IsInitialized() const = 0;

    virtual bool IsFinalized() const = 0;

  protected:
    EnvironmentManager() = default;
};

/// Holder for general data related to MPI (or suitable serial equivalents for non-MPI runs).
/** This class manages a registry of DataCommunicators, which can be used to perform MPI communication.
 *  @see DataCommunicator, MPIDataCommunicator.
 */
class KRATOS_API(KRATOS_CORE) ParallelEnvironment
{
  public:
    ///@name Type Definitions
    ///@{

    constexpr static bool MakeDefault = true;
    constexpr static bool DoNotMakeDefault = false;

    ///@}
    ///@name Access
    ///@{

    /// Retrieve a registered DataCommunicator instance.
    /** @param rName The name used to register the string. */
    static DataCommunicator& GetDataCommunicator(const std::string& rName);

    /// Retrieve the default DataCommunicator instance.
    static DataCommunicator& GetDefaultDataCommunicator();

    /// Set a new default DataCommunicator instance.
    /** @param rName The name the new default DataCommunicator was registered with. */
    static void SetDefaultDataCommunicator(const std::string& rName);

    /// Get the rank of the current process, as given by the default DataCommunicator.
    static int GetDefaultRank();

    /// Get the MPI Comm size, as given by the default DataCommunicator.
    static int GetDefaultSize();

    /// Get the MPI Comm size, as given by the default DataCommunicator.
    static std::string RetrieveRegisteredName(const DataCommunicator& rComm);

    ///@}
    ///@name Operations
    ///@{

    static void SetUpMPIEnvironment(EnvironmentManager::Pointer pEnvironmentManager);

    /**
     * @brief Registers the fill communicator factory
     * This method takes the provided fill communicator pointer factory and saves it to be used later on
     * @param CommunicatorFactory Factory function returning a pointer to a (parallel or serial) fill communicator
     */
    template<class TDataCommunicatorInputType>
    static void RegisterFillCommunicatorFactory(std::function<FillCommunicator::Pointer(ModelPart&, TDataCommunicatorInputType&)> FillCommunicatorFactory);

    /**
     * @brief Registers the communicator factory
     * This method takes the provided communicator pointer factory and saves it to be used later on
     * @param CommunicatorFactory Factory function returning a pointer to a (parrallel or serial) communicator
     */
    template<class TDataCommunicatorInputType>
    static void RegisterCommunicatorFactory(std::function<Communicator::UniquePointer(ModelPart&, TDataCommunicatorInputType&)> CommunicatorFactory);

    /**
     * @brief Create a fill communicator object
     * This method uses the previously registered fill communicator factory for the creation of a new fill communicator pointer
     * @param rModelPart Model part to which the fill communicator will be applied
     * @return FillCommunicator::Pointer Pointer to the new fill communicator instance
     */
    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please use \"CreateFillCommunicatorFromGlobalParallelism\" instead")
    static FillCommunicator::Pointer CreateFillCommunicator(ModelPart& rModelPart);

    /**
     * @brief Create a fill communicator object
     * This method uses the previously registered fill communicator factory for the creation of a new fill communicator pointer
     * @param rModelPart Model part to which the fill communicator will be applied
     * @param rDataCommunicatorName Name of the data communicator to be retrieved for the fill communicator construction
     * @return FillCommunicator::Pointer Pointer to the new fill communicator instance
     */
    static FillCommunicator::Pointer CreateFillCommunicatorFromGlobalParallelism(
      ModelPart& rModelPart,
      const std::string& rDataCommunicatorName);

    /**
     * @brief Create a fill communicator object
     * This method uses the previously registered fill communicator factory for the creation of a new fill communicator pointer
     * @param rModelPart Model part to which the fill communicator will be applied
     * @param rDataCommunicator reference to the data communicator to be used for the fill communicator construction
     * @return FillCommunicator::Pointer Pointer to the new fill communicator instance
     */
    static FillCommunicator::Pointer CreateFillCommunicatorFromGlobalParallelism(
      ModelPart& rModelPart,
      const DataCommunicator& rDataCommunicator);

    /**
     * @brief Create a Communicator object
     * This method uses the previously registered communicator factory for the creation of a new communicator pointer
     * @param rModelPart Model part required to retrieve the variables list from it
     * @param rDataCommunicatorName Name of the data communicator to be retrieved for the communicator construction
     * @return Communicator::UniquePointer Unique pointer to the new communicator
     */
    static Communicator::UniquePointer CreateCommunicatorFromGlobalParallelism(
      ModelPart& rModelPart,
      const std::string& rDataCommunicatorName);

    /**
     * @brief Create a Communicator object
     * This method uses the previously registered communicator factory for the creation of a new communicator pointer
     * @param rModelPart Model part required to retrieve the variables list from it
     * @param rDataCommunicator reference to the data communicator to be used for the communicator construction
     * @return Communicator::UniquePointer Unique pointer to the new communicator
     */
    static Communicator::UniquePointer CreateCommunicatorFromGlobalParallelism(
      ModelPart& rModelPart,
      const DataCommunicator& rDataCommunicator);

    /// Add a new DataCommunicator instance to the ParallelEnvironment.
    /** @param rName The name to be used to identify the DataCommunicator within ParallelEnvironment.
     *  @param rPrototype The DataCommunicator instance.
     *  @param Default If set to ParallelEnvironment::MakeDefault (a.k.a. true), the provided DataCommunicator will also be made default.
     */
    static void RegisterDataCommunicator(
        const std::string& rName,
        DataCommunicator::UniquePointer pPrototype,
        const bool Default = DoNotMakeDefault);

    /// Remove a DataCommunicator instance from the ParallelEnvironment.
    /** @param rName The name used to register the DataCommunicator within ParallelEnvironment.
     */
    static void UnregisterDataCommunicator(
        const std::string& rName);

    ///@}
    ///@name Inquiry
    ///@{

    /// Check if a DataCommunicator is registered as rName.
    static bool HasDataCommunicator(const std::string& rName);

    /// Get the registered name of the current default.
    /** This is a convenience function to help with temporarily changing the default DataCommunicator. */
    static std::string GetDefaultDataCommunicatorName();

    static bool MPIIsInitialized();

    static bool MPIIsFinalized();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    static std::string Info();

    /// Print information about this object.
    static void PrintInfo(std::ostream &rOStream);

    /// Print object's data.
    static void PrintData(std::ostream &rOStream);

    ///@}

  private:

    ///@name Private Life Cycle
    ///@{

    /// Internal constructor.
    ParallelEnvironment();

    /// Destructor
    ~ParallelEnvironment();

    ///@}
    ///@name Private Operations
    ///@{

    static void Create();

    void SetUpMPIEnvironmentDetail(EnvironmentManager::Pointer pEnvironmentManager);

    template<class TDataCommunicatorInputType>
    void RegisterFillCommunicatorFactoryDetail(std::function<FillCommunicator::Pointer(ModelPart&, TDataCommunicatorInputType&)> FillCommunicatorFactory);

    template<class TDataCommunicatorInputType>
    void RegisterCommunicatorFactoryDetail(std::function<Communicator::UniquePointer(ModelPart&, TDataCommunicatorInputType&)> CommunicatorFactory);

    void RegisterDataCommunicatorDetail(
        const std::string& Name,
        DataCommunicator::UniquePointer pPrototype,
        const bool Default = DoNotMakeDefault);

    void UnregisterDataCommunicatorDetail(const std::string& Name);

    DataCommunicator& GetDataCommunicatorDetail(const std::string& rName) const;

    DataCommunicator& GetDefaultDataCommunicatorDetail() const;

    void SetDefaultDataCommunicatorDetail(const std::string& rName);

    ///@}
    ///@name Private Inquiry
    ///@{

    bool HasDataCommunicatorDetail(const std::string& rName) const;

    bool MPIIsInitializedDetail() const;

    bool MPIIsFinalizedDetail() const;

    ///@}
    ///@name Private Access
    ///@{

    static ParallelEnvironment& GetInstance();

    void SetAsDefault(std::unordered_map<std::string, DataCommunicator::UniquePointer>::iterator& rThisCommunicator);

    ///@}
    ///@name Private Input and output
    ///@{

    /// Turn back information as a string.
    std::string InfoDetail() const;

    /// Print information about this object.
    void PrintInfoDetail(std::ostream &rOStream) const;

    /// Print object's data.
    void PrintDataDetail(std::ostream &rOStream) const;

    ///@}
    ///@name Member Variables
    ///@{

    std::unordered_map<std::string, DataCommunicator::UniquePointer> mDataCommunicators;

    std::unordered_map<std::string, DataCommunicator::UniquePointer>::iterator mDefaultCommunicator;

    std::function<FillCommunicator::Pointer(ModelPart&, const std::string&)> mFillCommunicatorStringFactory;

    std::function<FillCommunicator::Pointer(ModelPart&, const DataCommunicator&)> mFillCommunicatorReferenceFactory;

    std::function<Communicator::UniquePointer(ModelPart&, const std::string&)> mCommunicatorStringFactory;

    std::function<Communicator::UniquePointer(ModelPart&, const DataCommunicator&)> mCommunicatorReferenceFactory;

    int mDefaultRank;
    int mDefaultSize;

    EnvironmentManager::Pointer mpEnvironmentManager;

    static ParallelEnvironment* mpInstance;
    static bool mDestroyed;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ParallelEnvironment &operator=(ParallelEnvironment const &rOther) = delete;

    /// Copy constructor.
    ParallelEnvironment(ParallelEnvironment const &rOther) = delete;

    ///@}

}; // Class ParallelEnvironment

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_PARALLEL_ENVIRONMENT_H_INCLUDED  defined

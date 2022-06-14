# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from .wind_load_solver import WindLoadSolver
from KratosMultiphysics.CoSimulationApplication.utilities.data_communicator_utilities import GetRankZeroDataCommunicator

def Create(settings, model, solver_name):
    return SdofSolverWrapper(settings, model, solver_name)

class SdofSolverWrapper(CoSimulationSolverWrapper):
    """ This class implements a wrapper for an WindLoad solver to be used in CoSimulation
    """
    def __init__(self, settings, model, solver_name, model_part_name="Sdof"):
        super().__init__(settings, model, solver_name)

        self.mp = self.model.CreateModelPart(model_part_name)
        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 1

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        self._wind_solver = self._CreateWindLoadSolver(self.model, input_file_name)


    @classmethod
    def _CreateWindLoadSolver(cls, model, input_file_name):
        return WindLoadSolver(model, input_file_name)

    def Initialize(self):
        super().Initialize()
        self._wind_solver.Initialize()

    def AdvanceInTime(self, current_time):
        return self._wind_solver.AdvanceInTime(current_time)

    def SolveSolutionStep(self):

        self._wind_solver.SolveSolutionStep()

    def _GetDataCommunicator(self):
        # this solver does not support MPI
        return GetRankZeroDataCommunicator()

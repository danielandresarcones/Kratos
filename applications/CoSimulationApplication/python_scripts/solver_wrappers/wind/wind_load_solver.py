# CoSimulation imports
from asyncore import loop
import math
from random import random
from socketserver import ThreadingUnixStreamServer
from this import d
from time import time
import KratosMultiphysics
from KratosMultiphysics.CoSimulationApplication.function_callback_utility import GenericCallFunction
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Other imports
import numpy as np
import json
import os

class WindLoadSolver(MechanicalSolver):
    """ This class implements an WindLoad solver based on a MechanicalSolver
    It generates the loads produced by coherent wind at the nodes of a structure.
    """
    def __init__(self, model, input_name):

        # mimicing two constructors
        if isinstance(input_name, dict):
            parameters = KratosMultiphysics.Parameters(input_name)

        elif isinstance(input_name, str):
            if not input_name.endswith(".json"):
                input_name += ".json"

            with open(input_name,'r') as ProjectParameters:
                parameters = KratosMultiphysics.Parameters(ProjectParameters.read())

        else:
            raise Exception("The input has to be provided as a dict or a string")

        default_settings = self.GetDefaultParameters()

        RecursivelyValidateAndAssignDefaults(default_settings, parameters)

        self.model = model

        super().__init__(self.model, parameters)

        # Design parameters
        self.norm = parameters["design_parameters"]["norm"].GetString()
        self.peak_velocity = parameters["design_parameters"]["peak_velocity_kmh"].GetDouble()

        # TODO: Adapt for using the line height from model
        self.nominal_line_height = parameters["design_parameters"]["nominal_line_height"].GetDouble()
        self.reference_height = parameters["design_parameters"]["reference_height"].GetDouble()
        self.span_length = parameters["design_parameters"]["span_length"].GetDouble()
        self.terrain_rugosity_factor = parameters["design_parameters"]["terrain_rugosity_factor"].GetDouble()
        self.terraing_rugosity_equivalent_height = parameters["design_parameters"]["terraing_rugosity_equivalent_height"].GetDouble()
        self.decay_coefficient = parameters["design_parameters"]["decay_coefficient"].GetDouble()

        # Load parameters
        self.terrain_category = parameters["design_parameters"]["terrain_category"].GetInt()
        self.drag_coefficient = parameters["design_parameters"]["drag_coefficient"].GetDouble()
        self.conductor_diameter = parameters["design_parameters"]["conductor_diameter"].GetDouble()
        self.air_density = parameters["design_parameters"]["air_density"].GetDouble()

        # Wind model parameters
        self.turbulent_wind = parameters["wind_model_parameters"]["turbulent_wind"].GetBool()
        self.wind_model_name = parameters["wind_model_parameters"]["wind_model"].GetString()
        self.number_of_frequencies = parameters["wind_model_parameters"]["number_of_frequencies"].GetInt()
        self.frequency_herz = parameters["wind_model_parameters"]["frequency_in_herz"].GetBool()
        self.frequency_range_step = parameters["wind_model_parameters"]["frequency_range_step"].GetDouble()
        self.minimum_frequency = parameters["wind_model_parameters"]["minimum_frequency"].GetDouble()
        self.std_deviation_relation_resolution = parameters["wind_model_parameters"]["std_deviation_relation_resolution"].GetDouble()
        self.confidence_peak_velocity = parameters["wind_model_parameters"]["confidence_peak_velocity"].GetDouble()
        self.incidence_angle = parameters["wind_model_parameters"]["incidence_angle"].GetDouble()*math.pi/180.0
        self.ramp_up_time = parameters["wind_model_parameters"]["ramp_up_time"].GetDouble()

        self.wind_model = self.GetWindModelSND(self.wind_model_name, self.reference_height, self.nominal_line_height)

        # Output parameters
        self.write_output_file = parameters["output_parameters"]["write_output_file"].GetBool()
        self.file_name = parameters["output_parameters"]["file_name"].GetString()

        # Initialize geometry
        model_part_file = self.settings["model_import_settings"]["input_filename"].GetString()
        self.model_part_name = self.settings["model_part_name"].GetString()

        if self.model.HasModelPart(self.model_part_name):
            self.root_model_part = self.model.GetModelPart(self.model_part_name)
        else:
            self.root_model_part = self.model.CreateModelPart(self.model_part_name)
        
        self.AddVariables()
        self._add_dynamic_variables()
        self.root_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD)

        KratosMultiphysics.ModelPartIO(model_part_file).ReadModelPart(self.root_model_part)
        self.main_model_part = self.root_model_part.GetSubModelPart(self.model_part_name)

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type" : "wind_load_solver",
            "model_part_name" : "",
            "domain_size" : -1,
            "echo_level": 0,
            "buffer_size": 2,
            "analysis_type": "non_linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": ""
            },
            "design_parameters":{
                "norm": "euronorm",
                "peak_velocity_kmh": 140,
                "nominal_line_height": 30,
                "reference_height": 10,
                "terrain_rugosity_factor": 0.28,
                "terraing_rugosity_equivalent_height": 0.01,
                "decay_coefficient": 15,
                "span_length": 400,
                "terrain_category": 1,
                "drag_coefficient": 1.0,
                "conductor_diameter": 1E-2,
                "air_density": 1.225 
            },
            "wind_model_parameters":{
                "turbulent_wind": true,
                "wind_model": "kaimal",
                "number_of_frequencies": 400,
                "frequency_range_step": 5E-3,
                "minimum_frequency": 0.005,
                "std_deviation_relation_resolution": 10,
                "confidence_peak_velocity": 3,
                "frequency_in_herz": false,
                "incidence_angle": 90.0,
                "ramp_up_time": 0.0
            },
            "output_parameters":{
                "write_output_file": true,
                "file_name" : "output/results_wind.dat"
            },
            "time_stepping" : { },
            "volumetric_strain_dofs": false,
            "rotation_dofs": false,
            "pressure_dofs": false,
            "displacement_control": false,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "use_old_stiffness_in_first_iteration": false,
            "compute_reactions": true,
            "solving_strategy_settings": {
                "type" : "newton_raphson",
                "advanced_settings" : { }
            },
            "builder_and_solver_settings" : {
                "use_block_builder" : true,
                "use_lagrange_BS"   : false,
                "advanced_settings" : { }
            },
            "clear_storage": false,
            "move_mesh_flag": true,
            "multi_point_constraints_used": true,
            "convergence_criterion": "residual_criterion",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] = current_time

        return new_time



    def Initialize(self):

        # Initialize DOFs

        self.AddDofs()
        self._add_dynamic_dofs()

        # Initialize frequency space
        self.frequency_space = np.arange(0,self.number_of_frequencies)*self.frequency_range_step+self.minimum_frequency

        # Initialize distances matrix
        # TODO: This should be modified in each iteration, for now it is static. Also, nodes are supposed to be aligned,
        #       at the same heights y and z and generated with the generate_geometry.py tool
        # TODO: Uncouple the wind evaluation nodes from the geometry nodes

        self.distances_matrix = []
        self.number_of_windpoints = len(self.main_model_part.GetNodes())

        for first_node in self.main_model_part.GetNodes():
            distance_first_node = []
            for second_node in self.main_model_part.GetNodes():
                distance_first_node.append(abs(first_node.X-second_node.X))
            
            self.distances_matrix.append(distance_first_node)

        # Calculate design velocity
        self.CalculateObjectiveMeanVelocity()
        
        # Calculate PSDMatrix
        if self.turbulent_wind:
            self.PSD_Matrix = self.CaculatePSDMatrix(self.distances_matrix, self.frequency_space)

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
        pass

    def SolveSolutionStep(self):
        wind_velocity = self.CalculateVelocityVector()
        nodes_list = [node for node in self.main_model_part.GetNodes()]
        cable_velocity = [math.sqrt((node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)*math.sin(self.incidence_angle))**2 +
                          (node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)*math.cos(self.incidence_angle))**2) for node in nodes_list]
        relative_velocity = np.subtract(wind_velocity,cable_velocity)
        line_loads = [-self.CalculateLoad(velocity) for velocity in relative_velocity]
        node_loads = []
        for inode in range(len(nodes_list)):
            if inode == 0:
                node_distance = 0.5*abs(nodes_list[inode+1].X-nodes_list[inode].X)
            elif inode == len(nodes_list)-1:
                node_distance = 0.5*abs(nodes_list[inode].X-nodes_list[inode-1].X)
            else:
                node_distance = 0.5*abs(nodes_list[inode+1].X-nodes_list[inode-1].X)
                
            node_loads.append(line_loads[inode]*node_distance)



        for node, load in zip(nodes_list, node_loads):
            node.SetSolutionStepValue(KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD_Z, load * math.sin(self.incidence_angle))
            node.SetSolutionStepValue(KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD_X, load * math.cos(self.incidence_angle))

    def FinalizeSolutionStep(self):
        pass


    def CaculatePSDMatrix(self, distances_matrix, frequency_space):
        '''
        This function calculates the PSD Matrix ASSUMING SAME HEIGHT.
        '''

        #Initialize PSD matrix
        PSD_matrix = np.empty((len(distances_matrix),len(distances_matrix), self.number_of_frequencies))
        # TODO: Make it general for different heights

        # Calculate diagonal terms
        PSD_diagonal = self.CalculatePSDDiagonal(frequency_space, self.nominal_line_height)

        # Calculate off-diagonal terms
        for i in range(len(distances_matrix[:])):
            for j in range(i+1):
                if i==j:
                    PSD_matrix[i][j][:]=PSD_diagonal
                else:
                    PSD_matrix[i][j][:] = np.multiply(PSD_diagonal,np.exp(-(self.decay_coefficient*distances_matrix[i][j])/(4*math.pi*self.design_velocity)*np.absolute(frequency_space)))
                    PSD_matrix[j][i][:] = PSD_matrix[i][j][:]
        
        return PSD_matrix

    def CalculatePSDDiagonal(self, frequency_space, height):
        '''
        This function calculates the members of the diagonal of the PSD matrix for a discretized frequency space.
        '''
        S_diagonal = []
        deviation = self.design_velocity/np.log(height/self.terraing_rugosity_equivalent_height)
        for frequency in frequency_space:
            if not self.frequency_herz:
                frequency_herz = frequency/(2*math.pi)
                frequency_rad = frequency
            else:
                frequency_herz = frequency
                frequency_rad = frequency*2*math.pi

            S_diagonal.append(self.wind_model.CalculateSND(frequency_rad, self.design_velocity)*deviation**2/frequency_herz)

        return S_diagonal

    def CalculateVelocityVector(self):

        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if current_time > self.ramp_up_time:
            factor = 1.0
        else:
            factor = math.sin(math.pi/2.0*(current_time/self.ramp_up_time))
        
        if self.turbulent_wind:
            turbulent_velocity = self.CalculateTurbulentVelocityComponent(self.distances_matrix, current_time)

            return np.add([self.design_velocity*factor]*self.number_of_windpoints,turbulent_velocity)
        
        else:
            return [self.design_velocity*factor]*self.number_of_windpoints

    def CalculateTurbulentVelocityComponent(self, distances_matrix, time):

        u_turbulent = []
        PSD_matrix = self.PSD_Matrix
        H_matrices = [np.linalg.cholesky(PSD_matrix[:,:,ifrequency]) for ifrequency in range(self.number_of_frequencies)]

        for i in range(len(distances_matrix[:])):

            sum_1 = 0

            for j in range(i):
                sum_2 = 0
                fi = random()*2*math.pi

                for ifrequency in range(self.number_of_frequencies):

                    # Calculate H_matrix as the Cholesky decomposition of the PSD matrix
                    H_matrix = H_matrices[ifrequency]
                    h_ij = math.atan(np.imag(H_matrix[i,j]/np.real(H_matrix[i,j])))

                    sum_2 = sum_2 + abs(H_matrix[i,j])*math.cos(self.frequency_space[ifrequency]*time-h_ij+fi)
                
                sum_1 = sum_1 + sum_2

            u_turbulent.append(2*self.frequency_range_step*sum_1)
        
        return u_turbulent

    def CalculateObjectiveMeanVelocity(self):
        '''
            TODO: This function should estimate the design mean velocity for a given peak value according to a 
            distribution rule (e.g. 3-sigma rule). A sample of trial velocities is applied, the deviations are
            calculated and then an interpolation is performed to obtein the desired one.
            Set to 0.8 * peak velocity right now, to be changed in the future.
        '''

        self.design_velocity = 0.8 * self.peak_velocity


    def CalculateLoad(self, velocity):

        if self.norm == "euronorm":
            return self.CalculateLoadEuronorm(velocity)

        elif self.norm == "rlat":
            return self.CalculateLoadRLAT(velocity)

        else:
            raise Exception("Design according to norm {} not available. The norms implemented are euronorm and rlat.".format(self.norm))

    def CalculateLoadEuronorm(self, velocity):
        # TODO: Implement compensation factor for extreme winds according to Gumbel function for return times different to 50 years
        # TODO: Implement velocity design from different sources and wind maps
        # Using peak_velocity as design one

        # TODO: Recalculate air density for other parameters somewhere
        dynamic_load = 0.5*self.air_density*velocity**2

        gust_factor = (1.0+2.28/np.log(self.nominal_line_height/self.terraing_rugosity_equivalent_height))**2

        if self.terrain_category == 1:
            span_factor = 1.3-0.073*np.log(self.span_length)
        elif self.terrain_category == 2:
            span_factor = 1.3-0.082*np.log(self.span_length)
        elif self.terrain_category == 3:
            span_factor = 1.3-0.098*np.log(self.span_length)
        elif self.terrain_category == 4:
            span_factor = 1.3-0.110*np.log(self.span_length)
        else:
            raise Exception("Terrain categories must be 1-4.")
        
        force_conductor = dynamic_load * gust_factor * span_factor * self.drag_coefficient * self.conductor_diameter
        
        return force_conductor

    def CalculateLoadRLAT(self, velocity):

        if self.conductor_diameter < 1.6E-2:
            line_load = 60.0*(velocity/120.0)**2
        else:
            line_load = 50.0*(velocity/120.0)**2
        
        force_conductor = line_load * self.conductor_diameter         

        return force_conductor




    @staticmethod
    def GetWindModelSND(model, reference_height, height):
        wind_models = {
            "kaimal": KaimalModel,
            "davenport": DavenportModel,
            "von_karman": VonKarmanModel,
            "solari": SolariModel,
            "panofsky": PanofskyModel,
            "harris": HarrisModel
        }

        return wind_models[model](reference_height, height)

class WindModel(object):

    def __init__(self, reference_height, height):

        self.reference_height = reference_height
        self.height = height
        self.nominal_length = self.CalculateNominalLenght()

    def CalculateSND(self, frequency, velocity):
        self.frequency_L = self.CalculateFrequencyLength(frequency, velocity)

    def CalculateNominalLenght(self):
        return 0

    def CalculateFrequencyLength(self, frequency, velocity):

        return frequency * self.nominal_length / velocity         


class KaimalModel(WindModel):
    """ Kaimal spectral density model"""

    def CalculateNominalLenght(self):
        return 411 * (self.height/300)**(0.046+0.074*np.log(self.reference_height))

    def CalculateSND(self, frequency, velocity):
        super().CalculateSND(frequency, velocity)
        return 6.8*(self.frequency_L)/(1+10.2*self.frequency_L)**(5/3)

class DavenportModel(WindModel):
    """ Davenport spectral density model"""

    def CalculateNominalLenght(self):
        return 1200
    def CalculateSND(self, frequency, velocity):
        super().CalculateSND(frequency, velocity)
        return 2*(self.frequency_L)**2/(3*(1+self.frequency_L**2)**(5/3))

class VonKarmanModel(WindModel):
    """ Von Karman spectral density model"""

    def CalculateNominalLenght(self):
        return 300 * (self.height/300)**(0.046+0.074*np.log(self.reference_height))

    def CalculateSND(self, frequency, velocity):
        super().CalculateSND(frequency, velocity)
        return 4*(self.frequency_L)/(1+70.8*self.frequency_L**2)**(5/6)

class SolariModel(WindModel):
    """ Solari spectral density model"""

    def CalculateNominalLenght(self):
        return 411 * (self.height/300)**(0.046+0.074*np.log(self.reference_height))

    def CalculateSND(self, frequency, velocity):
        super().CalculateSND(frequency, velocity)
        return 6.868*(self.frequency_L)/(1+10.302*self.frequency_L)**(5/3)

class PanofskyModel(WindModel):
    """ Panofsky spectral density model"""

    def CalculateNominalLenght(self):
        return 411 * (self.height/300)**(0.046+0.074*np.log(self.reference_height))

    def CalculateSND(self, frequency, velocity):
        super().CalculateSND(frequency, velocity)
        return 1.1*(self.frequency_L)/(1+4*self.frequency_L)**(2)

class HarrisModel(WindModel):
    """ Harris spectral density model"""

    def CalculateNominalLenght(self):
        return 1800

    def CalculateSND(self, frequency, velocity):
        super().CalculateSND(frequency, velocity)
        return 6.8*(self.frequency_L)/(1+10.2*self.frequency_L)**(5/3)

def RecursivelyValidateAndAssignDefaults(defaults, settings):
    settings.ValidateAndAssignDefaults(defaults)

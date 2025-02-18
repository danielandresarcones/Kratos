import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class KinematicConstraintsTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "kinematic_constraints_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-3
        node = self.spheres_model_part.GetNode(1)
        velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
        angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
        if self.time > 0.19 and self.time < 0.2:
            expected_value = 0.0
            self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
            self.check_mark_1 = True
        elif self.time > 0.3095 and self.time < 0.31:
            expected_value = -1.0791000000000033
            self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
            self.check_mark_2 = True

        node = self.spheres_model_part.GetNode(2)
        velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
        if self.time > 0.25 and self.time < 0.3:
            expected_value = -10.0 * self.time
            self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
            self.check_mark_3 = True
        if self.time > 0.3095 and self.time < 0.31:
            expected_value = -0.09319499999999994
            self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
            self.check_mark_4 = True

        node = self.spheres_model_part.GetNode(3)
        velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
        angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
        if self.time < 0.1 and self.time > 0.05:
            expected_value = -5.0
            self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
            expected_value = 0.0
            self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
            expected_value = -10.0
            self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
            self.check_mark_5 = True

        node = self.spheres_model_part.GetNode(4)
        angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
        if self.time > 0.2495 and self.time < 0.25:
            expected_value = 0.21215864699077178
            self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
            self.check_mark_6 = True

        node = self.spheres_model_part.GetNode(5)
        velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
        angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
        if self.time > 0.2994 and self.time < 0.30:
            expected_value = 2.995000000000002
            self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
            expected_value = 2.990095000000002
            self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
            expected_value = 2.995000000000002
            self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
            expected_value = 2.995000000000002
            self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
            expected_value = 2.995000000000002
            self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
            expected_value = 2.995000000000002
            self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
            self.check_mark_7 = True

        node = self.rigid_face_model_part.GetNode(10)
        velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
        angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
        if self.time > 0.2995 and self.time < 0.30:
            expected_value = -0.464725
            self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
            expected_value = 1.0408537823387725
            self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
            expected_value = 9.191646120463801
            self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
            expected_value = 2.447989288623642
            self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
            expected_value = -0.6
            self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
            expected_value = 1.000
            self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
            self.check_mark_8 = True

        node = self.rigid_face_model_part.GetNode(20)
        velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
        angular_velocity = node.GetSolutionStepValue(Kratos.ANGULAR_VELOCITY)
        if self.time > 0.2995 and self.time < 0.30:
            expected_value = 2.4193720111091825
            self.CheckKinematicValues(velocity, 0, expected_value, tolerance)
            expected_value = -0.056663533294800894
            self.CheckKinematicValues(velocity, 1, expected_value, tolerance)
            expected_value = 0.6500061249830357
            self.CheckKinematicValues(velocity, 2, expected_value, tolerance)
            expected_value = 0.225
            self.CheckKinematicValues(angular_velocity, 0, expected_value, tolerance)
            expected_value = 1.1980000000000008
            self.CheckKinematicValues(angular_velocity, 1, expected_value, tolerance)
            expected_value = 2.695500000000002
            self.CheckKinematicValues(angular_velocity, 2, expected_value, tolerance)
            self.check_mark_9 = True

    def CheckKinematicValues(self, kinematic_variable, component, expected_value, tolerance):
        self.assertAlmostEqual(kinematic_variable[component], expected_value, delta=tolerance)

    def Finalize(self):
        self.assertTrue(self.check_mark_1)
        self.assertTrue(self.check_mark_2)
        self.assertTrue(self.check_mark_3)
        self.assertTrue(self.check_mark_4)
        self.assertTrue(self.check_mark_5)
        self.assertTrue(self.check_mark_6)
        self.assertTrue(self.check_mark_7)
        self.assertTrue(self.check_mark_8)
        self.assertTrue(self.check_mark_9)
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestKinematicConstraints(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def test_KinematicConstraints_1(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "kinematic_constraints_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(KinematicConstraintsTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()

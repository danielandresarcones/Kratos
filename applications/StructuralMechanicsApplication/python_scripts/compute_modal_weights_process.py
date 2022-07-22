import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return ComputeModalWeightsProcess(Model, settings["Parameters"])


class ComputeModalWeightsProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"           : "please_specify_model_part_name",
                "interval"                  : [0.0, 1e30],
                "write_output_file"    : true
            }
            """)

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString()):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30)
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.interval = KratosMultiphysics.Vector(2)
        self.interval[0] = settings["interval"][0].GetDouble()
        self.interval[1] = settings["interval"][1].GetDouble()
        self.write_output_file = settings["write_output_file"].GetBool()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):
                # Set drag output file name
                self.weights_filename = settings["model_part_name"].GetString() + "_modal_weigths.dat"

                # File creation to store the drag evolution
                with open(self.weights_filename, 'w') as file:
                    file.write(settings["model_part_name"].GetString() + " modal weights \n")
                    file.write("\n")
                    file.write("Time   Modal weights\n")
                    file.close()


    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and  (current_time < self.interval[1])):

            if (self.model_part.GetCommunicator().MyPID() == 0):

                if (self.write_output_file):
                    with open(self.weights_filename, 'a') as file:
                        output_str = str(current_time)
                        output_str += "   " + str(self.model_part.ProcessInfo.GetValue(KratosMultiphysics.StructuralMechanicsApplication.MODAL_WEIGHTS_VECTOR))
                        file.write(output_str)
                        file.close()

# PCMM: 

The PCMM Folder is meant to hold all PCMM projects

* This will be where we put our Ductal_Modeling project that we will turn into a github repo for publication

### Useful Tips

To compile scripts run the following commands in the PCMM directory or a specific project directory

    julia
    using PhysiCellModelManger
    initializeModelManager("Project_Name")

Next if you want to load an existing PhysiCell User Project: 

    importProject("../PhysiCell/user_projects/Duct_Dev"; dest=Dict("config" => "Duct_Dev_Updated", "custom_code" => "Duct_Dev_Updated", "rulesets_collection" => "Duct_Dev_Updated"))

To run scripts us the include command: 

    include("Duct_Dev/scripts/GenerateReport.jl")
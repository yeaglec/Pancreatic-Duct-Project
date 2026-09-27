## PhysiCell Ductal Project (Development)

This is my git-tracked PhysiCell project for developing my ductal model. It includes the following files: 

* Geometry.cpp - Code for membrane intialization and coordinate updates
* Membrane.cpp - Code for membrane related mechancis (strain, bending, etc.)
* Test.cpp - Test Suite for model mechanics and calulations
* custom.cpp - Staging File for PhysiCell wrapper and function hoopk-ups
* main.cpp - PhysiCell Compilation

PhysiCell_settings.xml - Model Specific Parameters

Some important paramters are the following: 
TODO 

### How to Use: 

In your PhysiCell Directory compile the ductal modeling project: 

    make clean
    make reset
    make load PROJ=Duct_Dev
    make 

Now in your Studio directory, run the following: 

    conda activate studio
    python bin/studio.py -c ../PhysiCell/user_projects/Duct_Dev/config/PhysiCell_settings.xml -e ../PhysiCell/project
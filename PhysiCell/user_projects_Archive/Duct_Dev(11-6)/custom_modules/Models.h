#ifndef __BM_models_h__
#define __BM_models_h__

#define _USE_MATH_DEFINES
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include <vector>
#include <string>
#include <tuple>

using namespace BioFVM; 
using namespace PhysiCell;


void cell_interactions_LSM(Cell* pCell,Phenotype& phenotype,double dt);
double distance_to_membrane(Cell* pCell,Phenotype& phenotype,double dt);
void test_enforce_boundary_repulsion(Cell* pCell,Phenotype& phenotype,double dt);
std::pair<double,double> basement_membrane_interactions_LSM(Cell* pCell);

void rebuild_signed_distance_field();
void update_basement_membrane_deformation(double dt);
void initialize_level_set_duct(std::vector<std::vector<double>> boundary_membrane_pts);




#endif
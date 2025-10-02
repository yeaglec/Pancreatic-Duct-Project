
#include "./Models.h" 
#include "./Utils.h"  
#include "./custom.h"
#include <cmath>
#include <cfloat>
#include <iostream>
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
// Models for the Level Set Method
// ___________________________________________________________________________________________________________________________

// ###################### Function that implements Cell-to-BM force ####################
void cell_interactions_LSM(Cell* pCell,Phenotype& phenotype, double dt){

    std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;
    double d = level_set_phi[i][j];
	

    double L = parameters.doubles("membrane_interaction_length"); // 500 now 
	double R = pCell->phenotype.geometry.radius;
	double de = d - (d < 0 ? -R : R);

    if(fabs(de) >= L) return;

	double cell_deadzone = parameters.doubles("cell_deadzone");
	
	if (fabs(de) <cell_deadzone) return;

	// Make this spring force

    auto grad   = level_set_gradient( i, j,level_set_phi, ls_dx, ls_dy);
    auto normal = level_set_normalize(grad);
    double sign     = (d < 0.0 ? +1.0 : -1.0);
    double strength = parameters.doubles("membrane_adhesion_strength"); //.001 right now
    double mag      = strength * fabs(de);

    pCell->velocity[0] += sign * mag * normal.first;
    pCell->velocity[1] += sign * mag * normal.second;
}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Level set functions for membrane interactions (voxel-based)
// ________________________________________________________________________________________________________________________

// Distance should just be the value of SDF at the voxel
double distance_to_membrane(Cell* pCell,Phenotype& phenotype,double dt){
    std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;

	return level_set_phi[i][j]; // SDF value at the voxel
}

void test_enforce_boundary_repulsion(Cell* pCell,Phenotype& phenotype,double dt){

	double cellx = pCell->position[0];
	double celly = pCell->position[1];

}

// ################ Function that computes the BM-to-Cell force ####################

std::pair<double,double> basement_membrane_interactions_LSM(Cell* pCell){       // TODO: Implement into update_basement_membrane_interactions{
	std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;

    double d = level_set_phi[i][j];
	// std::cout << "Distance to membrane: " << d << std::endl;
	double R = pCell->phenotype.geometry.radius;
	// std::cout << "Cell radius: " << R << std:: endl;
	double de = d - (d < 0 ? -R : R);                                    // de is the offseted distance to the BM
	// std::cout << "Effective distance to membrane: " << de << std::endl;

    double L = parameters.doubles("membrane_interaction_length");        // 500 now 
	double strength = parameters.doubles("membrane_spring_constant");  //.001 right now

	double BM_deadzone = parameters.doubles("membrane_deadzone");

	if (fabs(de) <BM_deadzone) return {0.0, 0.0}; 

    auto grad   = level_set_gradient( i, j,level_set_phi, ls_dx, ls_dy); // TThis is the gradient from the voxel center not the cell center
    auto normal = level_set_normalize(grad);
    double sign = (d < 0.0 ? +1.0 : -1.0);

    double mag = strength * fabs(de);    // ThisS is Hooke's law, F = kx

    double Fx = sign * mag * normal.first;
    double Fy = sign * mag * normal.second;

	return { Fx, Fy };
}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Functions for implementing deformations of basement membrane
// ________________________________________________________________________________________________________________________


void rebuild_signed_distance_field()
{
	// std::cout << "Rebuilding signed distance field..." << std::endl;
    // Mesh dimensions and spacing (already set in initialize duct code, just copied)
	
    int Nx = (int) level_set_phi.size();
    int Ny = (int) level_set_phi[0].size();

    // Build segment list from the current boundary points
    int Np = (int) boundary_membrane_pts.size();

    // For each grid cell, compute min distance to any segment
    for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
            double x = ls_xmin + (i + 0.5)*ls_dx;
			double y = ls_ymin + (j + 0.5)*ls_dy;
			auto [minDist, px, py, k, t] = project_point_onto_boundary(x, y);

            //Determine sign via is_inside() and write phi
            bool inside = is_inside(x, y, boundary_membrane_pts);
            level_set_phi[i][j] = inside ? -minDist : +minDist;
        }
    }
}

void update_basement_membrane_deformation(double dt){

	// std::cout << "Updating basement membrane deformation..." << std::endl;

	int Np = (int)boundary_membrane_pts.size();
	std::vector<std::pair<double,double>> node_forces(Np,{0,0});

	for (Cell* pCell : *all_cells){

		double cell_x = pCell->position[0];
		double cell_y = pCell->position[1];

		std::pair <double,double> force = basement_membrane_interactions_LSM(pCell);
		double Fx_cell = force.first;
		double Fy_cell = force.second;

		double Fx_BM = -Fx_cell;
		double Fy_BM = -Fy_cell;

		auto [best_dist, best_px, best_py, best_k_d, best_t] = project_point_onto_boundary(cell_x, cell_y);
		int best_k = static_cast<int>(std::round(best_k_d));

		// What this is doing:  take each cell’s force on the membrane, 
		//find which segment it hits, and split that tug between the two end‑nodes of that segment.

		int n1 = best_k, n2 = (best_k+1) % Np;
		double t_clamped = std::max(0.0, std::min(best_t, 1.0));  // Clamp t to [0,1]
		node_forces[n1].first += (1.0 - t_clamped) * Fx_BM;
		node_forces[n1].second += (1.0 - t_clamped) * Fy_BM;   // Split forces linearly based on t
		node_forces[n2].first += (t_clamped) * Fx_BM;
		node_forces[n2].second += (t_clamped) * Fy_BM;

	}

	for(int i=0; i<Np; i++) {
    boundary_membrane_pts[i][0] += node_forces[i].first  * dt;
    boundary_membrane_pts[i][1] += node_forces[i].second * dt;
	}

	// Rebuild the signed distance field after updating boundary points
	rebuild_signed_distance_field();

}

void initialize_level_set_duct(){
	std::cout << "Initializing level set function for the duct..." << std::endl;

	auto& mesh = get_default_microenvironment()->mesh;
	ls_xmin = mesh.bounding_box[0]; // bounding box is [xmin ymin zmin xmax ymax zmax]
	ls_ymin = mesh.bounding_box[1];
	ls_dx = mesh.dx;
	ls_dy = mesh.dy;     // Uniform grid spacing in x and y directions
	
	int Nx = (int)((mesh.bounding_box[3] - ls_xmin) / ls_dx);
	int Ny = (int)((mesh.bounding_box[4] - ls_ymin) / ls_dy); // Might need +1 here to include bb

	level_set_phi.assign(Nx, std::vector<double>(Ny, 0.0)); // Nx x Ny matrix initialized to zero

	// Creating BM shape
	int num_points = parameters.ints("membrane_num_points");
	double a = 300.0, b = 250.0;   // Length of ellipse axes
	double amp = 0.1;              // Amplitude of deformation
	int freq = 4;                 // Number of bumps

	double circle_radius = parameters.doubles("membrane_circle_radius");
	boundary_membrane_pts = generate_circle_boundary(circle_radius, num_points);

	int Np = (int)boundary_membrane_pts.size();
	initial_edge_length.resize(Np);

	// Store initial edge lengths for elasticity calculations
	for(int k = 0; k < Np; ++k)
	{
		int j = (k + 1) % Np;  
		double dx = boundary_membrane_pts[j][0] - boundary_membrane_pts[k][0];
		double dy = boundary_membrane_pts[j][1] - boundary_membrane_pts[k][1];
		initial_edge_length[k] = std::sqrt(dx*dx + dy*dy);
	}

	// High level overview: This double for loop gets the x and y coordinates of each voxel
	// Then computes the minimum distance to any segment (edge) of the shape

	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny; j++){

			// Coordinates at voxel center
			double x = ls_xmin + (i + 0.5)*ls_dx;
			double y = ls_ymin + (j + 0.5)*ls_dy;

			auto [minDist, px, py, k, t] = project_point_onto_boundary(x,y); // Project point onto boundary to get min distance
			bool inside = is_inside(x,y, boundary_membrane_pts);

			// Signed distance: negative inside, positive outside
			level_set_phi[i][j] = (inside ? -minDist : +minDist);
		}
	}
	std::cout << "Level set function initialized!!!" << std::endl;
}

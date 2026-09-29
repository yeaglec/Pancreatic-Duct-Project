#include "./Geometry.h"
#include "./Utils.h"
#include <cmath>
#include <cfloat>
#include <iostream>

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Membrane Architecture Helpers
// ________________________________________________________________________________________________________________________


/* 
Helper to test if point (x,y) is inside the polygon ##########
Determines if a point (x, y) is strictly inside a closed polygon using the 
ray-casting (even-odd) algorithm. 

Note on Edge case handling: for the intersection condition ((y1 > y) != (y2 > y))
Vertices: Strict inequalities ensure a horizontal ray passing exactly through 
a vertex is only counted if the connecting edge extends above the ray. This 
prevents double-counting when passing through a vertex.
Horizontal edges: Ignored safely as y1 and y2 will have the same relation to y.
*/
bool is_inside(double x, double y, const std::vector<std::vector<double>>& BM_pts) {

	// std::cout << "We are inside the is_inside function!!!" << std::endl;
	// std::cout << "Checking if point (" << x << ", " << y << ") is inside the polygon." << std::endl;
	
	bool inside = false;
	int Np = (int)BM_pts.size();

	for(int k=0; k<Np; ++k) {                                         

		double x1 = BM_pts[k][0], y1 = BM_pts[k][1];           
		double x2 = BM_pts[(k+1)%Np][0], y2 = BM_pts[(k+1)%Np][1];    

		// Check if horizontal ray intersects the edge between (x1,y1) and (x2,y2)
		if(((y1 > y) != (y2 > y)) ) {                                  
			double xint = x1 + (y - y1)*(x2 - x1)/(y2 - y1);           // Note here xint is the x-coordinate of the intersection point of the horizontal ray with the edge

			// Since the ray is cast in the +x direction, we only toggle the state 
            // if the intersection happens to the right of our test point.
			if(x < xint) inside = !inside;                             
		}
	}
return inside;
}

// ######### Function for finding to projection of a point onto the boundary #########
std::tuple<double, double, double, double, double> project_point_onto_boundary(double x, double y){   // TODO: Replace redundant projection code

	double best_dist = DBL_MAX;
	double best_px = 0.0, best_py = 0.0;
	int best_k = 0;
    double best_t = 0.0;
	int Np = boundary_membrane_pts.size();

	for (int k = 0; k < Np; ++k){

		int j = (k + 1) % Np;
		double x1 = boundary_membrane_pts[k][0], y1 = boundary_membrane_pts[k][1];
		double x2 = boundary_membrane_pts[j][0], y2 = boundary_membrane_pts[j][1];

		// Computing projection of vector c onto vector b

		double bx = x2 - x1, by = y2 - y1;                     // Vector b for segment
		double cx = x - x1, cy = y - y1;             // Vector c from first endpoint to cell position
		double b2 = bx*bx + by*by;
		double t = (b2 > 0.0 ? (cx*bx + cy*by) / b2 : 0.0);    // Just the projection formula: c \cdot b / ||b||^2

		double dist;
		double px, py;   // Our projection

		if (t <= 0.0) {
			// before first endpoint
			dist = sqrt(cx*cx + cy*cy);
			px = x1; py = y1;
		} else if (t >= 1.0) {
			// after second endpoint
			double ux = x - x2, uy = y - y2;
			dist = sqrt(ux*ux + uy*uy);
			px = x2; py = y2;
		} else {
			// projection between endpoints
			px = x1 + t * bx;
			py = y1 + t * by;
			double ux = x - px, uy = y - py;
			dist = sqrt(ux*ux + uy*uy);
		}

		if (dist < best_dist) {
			best_dist = dist;
			best_k = k;
			best_t = t;
			best_px = px;
			best_py = py;
		}
	}

	return {best_dist, best_px, best_py, best_k, best_t};

}

// Rudimentary function to generate duct shape
std::vector<std::vector<double>> generate_boundary_shape(double a, double b, double amp, int freq){

	int num_points = parameters.ints("membrane_num_points");
	std::vector<std::vector<double>> pts;
	pts.reserve(num_points);
    
	// Generating points for a deformed ellipse like shape
    for (int i = 0; i < num_points; i++) {
        
        double theta = 2.0 * M_PI * i / (num_points);
        
        double r_x = a * (1.0 + amp * cos(freq * theta));
        double r_y = b * (1.0 + amp * sin(freq * theta));
        double x = r_x * cos(theta);
        double y = r_y * sin(theta);

        pts.push_back({ x, y, 0.0 });
    }
    return pts;
}

// Cancer cells are placed around cell 0 in the order 0, 1, -1, 2, -2, ...
bool is_cancer_cell(int i, int num_cells, int num_cancer){
	int d = std::min(i, num_cells - i);                     // ring distance from cell 0
	int rank = (d == 0) ? 0 : (i == d ? 2*d - 1 : 2*d);     // placement order
	return rank < num_cancer;
}

// Generating points for a deformed ellipse like shape
void generate_boundary_cells(double a, double b, double amp, int freq, std::string type, double dis, int num_cells, int num_cancer){
	
	Cell_Definition* pBM_def = cell_definitions_by_name[type];   //cell_definitions_by_name[ type_name ] = pCD; 
	

    for (int i = 0; i < num_cells; i++) {
        
        double theta = 2.0 * M_PI * i / (num_cells);
        
        double r_x = a * (1.0 + amp * cos(freq * theta));
        double r_y = b * (1.0 + amp * sin(freq * theta));
        double x = r_x * cos(theta);
        double y = r_y * sin(theta);

		// compute radial distance and unit‐radial direction
		double r_norm = std::sqrt(x*x + y*y);
		double nx = x / r_norm;    // outward radial unit vector
		double ny = y / r_norm;

		// step back along the normal by ep_dis
		double xi = x - dis * nx;
		double yi = y - dis * ny;

		// NOTE: cancer cells currently use the "CAF" cell definition (naming to be fixed)
		Cell* pC = nullptr;
		if (is_cancer_cell(i, num_cells, num_cancer)) pC = create_cell( *cell_definitions_by_name["CAF"] );
		else pC = create_cell( *pBM_def );
		
		if (parameters.ints("number_total_cells") == 1) pC->assign_position( { parameters.doubles("x"), parameters.doubles("y"), 0.0 } );
		else pC->assign_position( {xi, yi, 0.0 } );
		
		if(i==0)pC->phenotype.cycle.data.exit_rate(0) = parameters.doubles("proliferation_exit_rate");

	
	}
}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Intialization Functions
// ________________________________________________________________________________________________________________________

std::vector<std::vector<double>> generate_circle_boundary(){
	double radius = parameters.doubles("membrane_circle_radius");
	int num_points = parameters.ints("membrane_num_points");
	std::vector<std::vector<double>> pts;	
	pts.reserve(num_points);

	for (int i = 0; i < num_points; i++) {
		double theta = 2.0 * M_PI * i / num_points;
		double x = radius * cos(theta);
		double y = radius * sin(theta);
		pts.push_back({ x, y, 0.0 });
	}
	return pts;
}

void generate_circle_cells(){
	
		int num_ep = parameters.ints("number_total_cells");
		int num_cancer = parameters.ints("num_cancer_cells");
		Cell_Definition* pBM_def = cell_definitions_by_index[0];
		double ep_dis = parameters.doubles("ep_displacement");

		for (int i=0; i<num_ep; i++) {
			double theta = 2.0 * M_PI * i / num_ep;
			double radius = parameters.doubles("membrane_circle_radius");
			double x = radius * cos(theta);
			double y = radius * sin(theta);

			// compute radial distance and unit‐radial direction
			double r_norm = std::sqrt(x*x + y*y);
			double nx = x / r_norm;    // outward radial unit vector
			double ny = y / r_norm;

			// step back along the normal by ep_dis
			double xi = x - ep_dis * nx;
			double yi = y - ep_dis * ny;

			// NOTE: cancer cells currently use the "CAF" cell definition (naming to be fixed)
			Cell* pC = nullptr;
			if (is_cancer_cell(i, num_ep, num_cancer)) pC = create_cell( *cell_definitions_by_name["CAF"] );
			else pC = create_cell( *pBM_def );

			// Turn on proliferation for first cell

			if( parameters.ints("number_total_cells") > 1 ){
				pC->assign_position( { xi, yi, 0.0 } );
			}
			else{
				pC->assign_position( { parameters.doubles("x"), parameters.doubles("y"), 0.0 } );
			}
			
			if (i==0){
				std::cout << "Setting first cell to proliferate" << std::endl;
				pC->phenotype.cycle.data.exit_rate(0) = parameters.doubles("proliferation_exit_rate");
			}

		}
	}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Initial Membrane Setups (toggled by user parameter "membrane_shape")
// ________________________________________________________________________________________________________________________

// 0: Deformed ellipse (default)
void setup_membrane_default(){
	double a = 300.0, b = 250.0;
	double amp = 0.1;              // Amplitude of deformation
	int freq = 4;
	int num_ep = parameters.ints("number_total_cells");
	int num_cancer = parameters.ints("num_cancer_cells");
	int num_caf = parameters.ints("number_CAF_cells");      // leftover: actual CAFs placed just outside the membrane
	double ep_dis = parameters.doubles("ep_displacement");

	boundary_membrane_pts = generate_boundary_shape(a, b, amp, freq);
	generate_boundary_cells(a, b, amp, freq, "Epithelial", ep_dis, num_ep, num_cancer);
	generate_boundary_cells(a, b, amp, freq, "CAF", -5, num_caf, 0);
}

// 1: Star (hard coded example)
void setup_membrane_star(){
	double a = 250.0, b = 250.0;
	double amp = 0.2;              // Amplitude of deformation
	int freq = 5;
	int num_ep = parameters.ints("number_total_cells");
	int num_cancer = parameters.ints("num_cancer_cells");
	int num_caf = parameters.ints("number_CAF_cells");      // leftover: actual CAFs placed just outside the membrane
	double ep_dis = parameters.doubles("ep_displacement");

	boundary_membrane_pts = generate_boundary_shape(a, b, amp, freq);
	generate_boundary_cells(a, b, amp, freq, "Epithelial", ep_dis, num_ep, num_cancer);
	generate_boundary_cells(a, b, amp, freq, "CAF", -5, num_caf, 0);
}

// 2: Circle (radius from "membrane_circle_radius")
void setup_membrane_circle(){
	boundary_membrane_pts = generate_circle_boundary();
	generate_circle_cells();
}

// TODO CLEAN: MOVE TO UTILS
void clustered_cell(double a, double b, double amp, int freq, int num_points, double center_x, double center_y, double radius, std::string type)
{
	Cell_Definition* pBM_def = cell_definitions_by_name[type];   //cell_definitions_by_name[ type_name ] = pCD; 
	std::cout << "Creating clustered BM cells..." << std::endl;
	std::vector<std::vector<double>> pts;
	pts.reserve(num_points);


    for (int j = 0; j < num_points; j++) {
        // Place cells randomly within a disc for a more real look
        double r = radius * sqrt(UniformRandom()); // sqrt for uniform area distribution
        double theta = 2.0 * M_PI * UniformRandom();
        
        double cx = center_x + r * cos(theta);
        double cy = center_y + r * sin(theta);
        
        Cell* pC2 = create_cell(*pBM_def);
        pC2->assign_position({cx, cy, 0.0}); 
    }

}
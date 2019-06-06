/*
	Created by Yuhan Ping from HKU at 2019-06-04.
	The project is about tooth completion.
	Under development.
	Lateset updated on 2019-06-06.
*/


#include <iostream>

#include "Header/io.h"
#include "Header/boundary_construction.h"

#include <igl/fit_plane.h>


using namespace Eigen;
using namespace std;


// ---------------pre-define parameters------------------
MatrixXd V_plane, V_total, V1_boundary_L, V1_boundary_R, V2_boundary_L, V2_boundary_R;

RowVectorXd L1, L2;
RowVector3d N, C; 


string first_mesh = "F:/StudyMaterials/HKU/RA/libigl-example-project/32.off";
string second_mesh = "F:/StudyMaterials/HKU/RA/libigl-example-project/33.off";
/*
	N: 1 * 3 vector, the normal of fitted plane
	C: 1 * 3 vector, a point lies in the plane
*/
//--------------------------------------------------------

int main(int argc, char *argv[])
{
	// load the mesh, currently load one mesh for test, load two mesh in the future for dealing with neighbor intersection
	load_mesh(first_mesh, V1, F1);

	// identify the boundary
	calc_average_edge_length(V1, F1);

	// create interactive straight line

	// TODO: triangulate the hole 

	// TODO: refinement

	// TODO: fairing with the surface conditions

	// Launch the viewer
	viewer_display(V1, F1);

	return 0;
}




/*
	Created by Yuhan Ping from HKU at 2019-06-04.
	The project is about tooth completion.
	Under development.
	Lateset updated on 2019-06-04.
*/


#include <iostream>
#include <igl/readOFF.h>
#include "Header/io.h"
#include <igl/boundary_loop.h>
#include <igl/fit_plane.h>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/ml.hpp>

using namespace Eigen;
using namespace std;
using namespace cv;
using namespace cv::ml;

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
	//igl::boundary_loop(F1, L1);
	//cout << "boundary vertices idx " << L1 << endl;

	// create interactive straight line

	// TODO: triangulate the hole 

	// TODO: refinement

	// TODO: fairing with the surface conditions

	// Launch the viewer
	viewer_display(V1, F1);

	return 0;
}




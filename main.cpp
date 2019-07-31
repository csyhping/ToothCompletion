/*
	Created by Yuhan Ping from HKU at 2019-06-04.
	The project is about tooth completion.
	Under development.
	Lateset updated on 2019-06-11.
*/


#include <iostream>

#include "Header/io.h"
#include "Header/boundary_construction.h"
#include "Header/plane.h"




using namespace Eigen;
using namespace std;



//--------------------------------------------------------

int main(int argc, char *argv[])
{
	if (argc != 2) {
		cout << "[Usage] ToothCompletion [path]" << endl;
		exit(-1);
	}

	// setup path
	inputmesh = argv[1];
	int nPos = inputmesh.find(".off");
	if (nPos != -1) {
		prefair_R = inputmesh.substr(0, nPos) + "\prefair_R.dmat";
		prefair_L = inputmesh.substr(0, nPos) + "\prefair_L.dmat";
		postfair_X = inputmesh.substr(0, nPos) + "\postfair_X.dmat";
		prefair_file = inputmesh.substr(0, nPos) + "\prefair.off";
		postfair_file = inputmesh.substr(0, nPos) + "\postfair.off";

	}


	// load the mesh, currently load one mesh for test, load two mesh in the future for dealing with neighbor intersection
 	load_mesh(inputmesh, V1, F1);

	// identify the boundary
	calc_average_edge_length(V1, F1);

	// visualize the boundary

	// create interactive straight line

	// TODO: triangulate the hole 

	// TODO: refinement

	// TODO: fairing with the surface conditions

	// Launch the viewer
	viewer_display(V1, F1);

	return 0;
}




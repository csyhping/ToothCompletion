#include "Header/plane.h"


void get_plane(Eigen::MatrixXd &Hole_vertex, Eigen::RowVector3d &N, Eigen::RowVector3d &C) {

	igl::fit_plane(Hole_vertex, N, C);
	std::cout << "the normal of the plane [" << N <<"]"<< std::endl;
	std::cout << "one point on the plane [" << C << "]" << std::endl;
}
void project_hole_vertex_to_plane();
void project_hole_vertex_back();
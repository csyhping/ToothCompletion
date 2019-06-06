#include "Header/boundary_construction.h"

Eigen::MatrixXd Color_per_vertex;
double avg_length; // store the average edge length

double calc_average_edge_length(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

	avg_length = igl::avg_edge_length(V, F);
	std::cout << "average edge length " << avg_length << std::endl;
	return avg_length;
}

void create_vertex_on_line(Eigen::RowVector3d &select_v1, Eigen::RowVector3d &select_v2) {

}
/*
	// visualization test, color matrix initialization
	Color_per_vertex = Eigen::MatrixXd::Constant(V.rows(), 3, 1);

	Eigen::RowVectorXd Boundary_loop; // store the ordered boundary loop of vertex idx, [v0, v1, v5, v6,...]
	igl::boundary_loop(F, Boundary_loop);
	std::cout << "boundary loop " << Boundary_loop << std::endl;
	std::cout << "boundary row col " << Boundary_loop.rows() << " " << Boundary_loop.cols() << std::endl;
	// modify color according for boundary vertex
	for (int i = 0; i < Boundary_loop.cols(); i++) {
		std::cout << "change color at " << Boundary_loop(0, i) << std::endl;
		Color_per_vertex.row(Boundary_loop(0, i)) << 1, 0, 0;
	}
*/
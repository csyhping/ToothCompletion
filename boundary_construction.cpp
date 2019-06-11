#include "Header/boundary_construction.h"

Eigen::MatrixXd Color_per_vertex;
Eigen::MatrixXd Hole_vertex;
double avg_length; // store the average edge length
int count; // how many new vertex

double calc_average_edge_length(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

	avg_length = igl::avg_edge_length(V, F);
	std::cout << "average edge length " << avg_length << std::endl;
	return avg_length;
}

void create_vertex_on_line(Eigen::RowVector3d &select_v1, Eigen::RowVector3d &select_v2, Eigen::MatrixXd &New_v_on_line) {

	// calculate line length
	Eigen::MatrixXd length;
	Eigen::MatrixXd v(2, 3);
	Eigen::MatrixXi f = (Eigen::MatrixXi(1, 2) <<
		0, 1).finished();
	v.row(0) = select_v1;
	v.row(1) = select_v2;
	igl::edge_lengths(v, f, length);
	std::cout << "length = " << length << std::endl;

	// how many new vertices
	count = length(0, 0) / avg_length;
	std::cout << "#new vertices to add " << count << std::endl;

	// calc coordinates of these vertices
	New_v_on_line.resize(count, 3); // resize according to #new_vertex

	for (int i = 0; i < count; i++) {
		New_v_on_line(i, 0) = select_v1(0, 0) + (-1)*(i + 1)*((select_v1(0, 0) - select_v2(0, 0)) / (count + 1));
		New_v_on_line(i, 1) = select_v1(0, 1) + (-1)*(i + 1)*((select_v1(0, 1) - select_v2(0, 1)) / (count + 1));
		New_v_on_line(i, 2) = select_v1(0, 2) + (-1)*(i + 1)*((select_v1(0, 2) - select_v2(0, 2)) / (count + 1));
	}
	std::cout << "new vertices coordinates " << New_v_on_line << std::endl;
	
}

void get_hole_boundary(Eigen::MatrixXd &color_per_v, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::RowVector3d &select_v1, Eigen::RowVector3d &select_v2, Eigen::MatrixXd &New_v_on_line) {
	
	//// visualization test, color matrix initialization
	//Color_per_vertex = Eigen::MatrixXd::Constant(V.rows(), 3, 1);

	Eigen::RowVectorXd Boundary_loop; // store the ordered boundary loop of vertex idx, [v0, v1, v5, v6,...]
	igl::boundary_loop(F, Boundary_loop);
	std::cout << "boundary loop " << Boundary_loop << std::endl;
	std::cout << "boundary row col " << Boundary_loop.rows() << " " << Boundary_loop.cols() << std::endl;
	//
	//// modify color according for boundary vertex to visualize original boundary 
	//for (int i = 0; i < Boundary_loop.cols(); i++) {
	//	std::cout << "change color at " << Boundary_loop(0, i) << std::endl;
	//	Color_per_vertex.row(Boundary_loop(0, i)) << 1, 0, 0;
	//}

	// construct new hole boundary 
	int size = 1 + 1 + count;
	int orignal_boundary_valid = 0;
	std::cout << "size = " << size << std::endl;
	double mid_x_of_select_v = 0.5*(select_v1(0, 0) + select_v2(0, 0));
	std::cout << "mid x = " << mid_x_of_select_v << std::endl;
	for (int i = 0; i < Boundary_loop.cols(); i++) {
		if (V(Boundary_loop(0, i), 0) <= mid_x_of_select_v) {
			std::cout << "v [" <<Boundary_loop(0,i)<<"]  is closer to v1 "<< std::endl;
			//// if close to select_v1
			//if (V(Boundary_loop(0, i), 2) > select_v1(0, 2)) {
			//	std::cout << "an original v is valid idx = " << Boundary_loop(i, 0) << std::endl;
			//	std::cout << "the coordinates of v is [" << V(Boundary_loop(0, i)) << std::endl;

			//	orignal_boundary_valid += 1;
			//	// if z > z_select_v1, add to hole_vertex
			//	Hole_vertex.resize(orignal_boundary_valid, 3);
			//	Hole_vertex.row(orignal_boundary_valid - 1) = V.row(Boundary_loop(0, i));
			//	std::cout << "hole vertex update " << Hole_vertex << std::endl;
			//	getchar();
			//}
		}
		else {
			// if close to select v2
			std::cout << "v [" << Boundary_loop(0, i) << "]  is closer to v2 " << std::endl;

		}
	}
}

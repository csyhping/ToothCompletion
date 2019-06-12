#include "Header/boundary_construction.h"

Eigen::MatrixXd Color_per_vertex;
Eigen::MatrixXd Hole_vertex_R, Hole_vertex_L;
double avg_length; // store the average edge length


double calc_average_edge_length(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {

	avg_length = igl::avg_edge_length(V, F);
	std::cout << "average edge length " << avg_length << std::endl;
	return avg_length;
}

void create_vertex_on_line(Eigen::RowVector3d &select_v1, Eigen::RowVector3d &select_v2, Eigen::MatrixXd &New_v_on_line, int &count) {

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

void get_hole_boundary(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::RowVector3d &select_v1, Eigen::RowVector3d &select_v2, Eigen::MatrixXd &New_v_on_line, int &idx_v1, int &idx_v2, int &count) {
	
	// visualization test, color matrix initialization
	Color_per_vertex = Eigen::MatrixXd::Constant(V.rows(), 3, 1);

	Eigen::RowVectorXd Boundary_loop; // store the ordered boundary loop of vertex idx, [v0, v1, v5, v6,...]
	igl::boundary_loop(F, Boundary_loop);
	std::cout << "boundary loop " << Boundary_loop << std::endl;
	std::cout << "boundary row col " << Boundary_loop.rows() << " " << Boundary_loop.cols() << std::endl;
	
	// construct a look up matrix to store the index of each element
	Eigen::RowVectorXi Pos_boundary;
	int max_v_idx = Boundary_loop.maxCoeff();
	std::cout << "max v index " << max_v_idx << std::endl;
	Pos_boundary.resize(max_v_idx + 1);
	Pos_boundary.setConstant(-1);
	std::cout << "pos boundary size " << Pos_boundary.rows() << " " << Pos_boundary.cols() << std::endl;
	for (int i = 0; i < Boundary_loop.cols(); i++) {
		Pos_boundary(0, Boundary_loop(0, i)) = i;
		std::cout << "v_" << Boundary_loop(0, i) << " at " << Pos_boundary(0, Boundary_loop(0, i)) << std::endl;
	}

	
	// pick a subsequent of vertices from boundary_loop from select_v1 to select_v2

	// construct new hole boundary 
	int size_boundary_R = Pos_boundary(0, idx_v2) - Pos_boundary(0, idx_v1) + 1 + count;
	std::cout << "size boundary R " << size_boundary_R << std::endl;
	Hole_vertex_R.resize(size_boundary_R, 3);
	std::cout << "hole boundary v rows cols " << Hole_vertex_R.rows() << " " << Hole_vertex_R.cols() << std::endl;

	int pos_boundary_loop_count = 0; // counter in boundary_loop
	int pos_new_v_on_line_count = 0; // counter in new_v_on_line

	for (int i = 0; i < size_boundary_R; i++) {
		if (pos_boundary_loop_count > (Pos_boundary(0, idx_v2) - Pos_boundary(0, idx_v1))) {
			// new created vertex on line 
			Hole_vertex_R.row(i) = New_v_on_line.row(pos_new_v_on_line_count);
			//std::cout << "load new v_" << pos_new_v_on_line_count << std::endl;
			pos_new_v_on_line_count += 1;
		}
		else
		{
			// original boundary vertices
			Hole_vertex_R.row(i) = V.row(Boundary_loop(0, Pos_boundary(0, idx_v1) + pos_boundary_loop_count));
			Color_per_vertex.row(Boundary_loop(0, Pos_boundary(0, idx_v1) + pos_boundary_loop_count)) << 1, 0, 0;
			std::cout << "load v_" << Boundary_loop(0, Pos_boundary(0, idx_v1) + pos_boundary_loop_count) << std::endl;
			pos_boundary_loop_count += 1;
		}

	}
}

#include "Header/boundary_construction.h"
#include "Header/io.h"

Eigen::MatrixXd Color_per_vertex;
Eigen::RowVectorXi Boundary_loop;
Eigen::RowVectorXi Pos_boundary;
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

void get_pos_boundary(Eigen::MatrixXi &F) {
	
	// get the boundary loop
	igl::boundary_loop(F, Boundary_loop);
	std::cout << "boundary loop " << Boundary_loop << std::endl;
	std::cout << "Totally boundary v " << Boundary_loop.cols() << std::endl;

	// construct a look up matrix to store the index of each element
	int max_v_idx = Boundary_loop.maxCoeff();
	std::cout << "max v idx = " << max_v_idx << std::endl;
	Pos_boundary.resize(max_v_idx+1);
	Pos_boundary.setConstant(-1);
	std::cout << "pos boundary size " << Pos_boundary.rows() << " " << Pos_boundary.cols() << std::endl;
	for (int i = 0; i < Boundary_loop.cols(); i++) {
		Pos_boundary(Boundary_loop(i)) = i;
		//std::cout << "v_" << Boundary_loop(i) << " at " << Pos_boundary(Boundary_loop(i)) << std::endl;
	}
}

void get_hole_boundary(Eigen::MatrixXd &V, Eigen::MatrixXd &New_v_on_line, int &idx_v1, int &idx_v2, int &count, Eigen::MatrixXd &Hole_vertex, Eigen::RowVectorXi &Hole_vertex_idx) {
	std::cout << "idx v1 , idx v2 =  " << idx_v1 << " , " << idx_v2 << std::endl;
	// construct new hole boundary 
	if (Pos_boundary(0, idx_v1) < Pos_boundary(0, idx_v2)) {
		// construct hole boundary if idx_v1 and idx_v2 are continues(no matter left or right side, continuous means they do not cover the start point of the boundary loop)
		int size_boundary_R;
		size_boundary_R = Pos_boundary(0, idx_v2) - Pos_boundary(0, idx_v1) + 1 + count;
		std::cout << "[INFO] the two selected vertex idx not cover origin idx " << std::endl;
		std::cout << "size boundary R " << size_boundary_R << std::endl;
		Hole_vertex.resize(size_boundary_R, 3);
		std::cout << "hole boundary v rows cols " << Hole_vertex.rows() << " " << Hole_vertex.cols() << std::endl;

		Hole_vertex_idx.resize(Pos_boundary(0, idx_v2) - Pos_boundary(0, idx_v1) + 1);
		std::cout << "hole vertex idx rows cols " << Hole_vertex_idx.rows() << " " << Hole_vertex_idx.cols() << std::endl;

		int pos_boundary_loop_count = 0; // counter in boundary_loop
		int pos_new_v_on_line_count = 0; // counter in new_v_on_line

		for (int i = 0; i < size_boundary_R; i++) {
			if (pos_boundary_loop_count > (Pos_boundary(0, idx_v2) - Pos_boundary(0, idx_v1))) {
				// new created vertex on line 
				Hole_vertex.row(i) = New_v_on_line.row(pos_new_v_on_line_count);
				//std::cout << "load new v_" << pos_new_v_on_line_count << std::endl;
				pos_new_v_on_line_count += 1;
			}
			else
			{
				// original boundary vertices
				Hole_vertex.row(i) = V.row(Boundary_loop(0, Pos_boundary(0, idx_v1) + pos_boundary_loop_count));
				Hole_vertex_idx(0,i) = Boundary_loop(0, Pos_boundary(0, idx_v1) + pos_boundary_loop_count);
				Color_per_vertex.row(Boundary_loop(0, Pos_boundary(0, idx_v1) + pos_boundary_loop_count)) << 1, 0, 0;
				pos_boundary_loop_count += 1;
			}

		}

	}
	else {
		// construct hole boundary if idx_v1 and idx_v2 covers the start point of the boundary loop, this only happens for left side
		int size_boundary_L;
		size_boundary_L = Boundary_loop.cols() +Pos_boundary(idx_v2)-Pos_boundary(idx_v1) + 1 + count;
		std::cout << "[INFO] the two selected vertex idx cover origin idx " << std::endl;
		cover_origin = 1;
		count_cover_new = count;
		std::cout << "size boundary L " << size_boundary_L << std::endl;
		Hole_vertex.resize(size_boundary_L, 3);
		std::cout << "hole boundary v rows cols " << Hole_vertex.rows() << " " << Hole_vertex.cols() << std::endl;
		Hole_vertex_idx.resize(Boundary_loop.cols() + Pos_boundary(idx_v2) - Pos_boundary(idx_v1) + 1);
		std::cout << "hole vertex idx rows cols " << Hole_vertex_idx.rows() << " " << Hole_vertex_idx.cols() << std::endl;

		int pos_new_v_on_line_count = 0; // counter in new_v_on_line
		int pos_v3_to_end = 0; // counter in boundary_loop from v3 to end point
		int pos_hole_idx = 0;// counter in hole_vertex_idx

		for (int i = 0; i < size_boundary_L; i++) {
			if (i <= Pos_boundary(idx_v2)) {
				// from start point to select_v4
				Hole_vertex.row(i) = V.row(Boundary_loop(i));
				Hole_vertex_idx(i) = Boundary_loop(i);
				pos_hole_idx += 1;
				count_cover_hole_part_left += 1;
				Color_per_vertex.row(Boundary_loop(i)) << 1, 0, 0;
				//pos_boundary_loop_count += 1;
			}
			else if (i > Pos_boundary(idx_v2) && pos_new_v_on_line_count < count) {
				// new created vertices
				Hole_vertex.row(i) = New_v_on_line.row(pos_new_v_on_line_count);
				pos_new_v_on_line_count += 1;
			}
			else if (i > Pos_boundary(idx_v2) && pos_new_v_on_line_count == count) {
				// from select_v3 to start point
				count_cover_hole_part_right += 1;
				Hole_vertex.row(i) = V.row(Boundary_loop(Pos_boundary(idx_v1) + pos_v3_to_end));
				Hole_vertex_idx(pos_hole_idx) = Boundary_loop(Pos_boundary(idx_v1) + pos_v3_to_end);
				Color_per_vertex.row(Boundary_loop(Pos_boundary(idx_v1) + pos_v3_to_end)) << 1, 0, 0;
				pos_v3_to_end += 1;
			}
		}

	}
	//std::cout << Hole_vertex << std::endl;
	std::cout << "hole vertex idx " << Hole_vertex_idx << std::endl;

}



#include "Header/plane.h"


void get_plane(Eigen::MatrixXd &Hole_vertex, Eigen::RowVector3d &N, Eigen::RowVector3d &C) {

	igl::fit_plane(Hole_vertex, N, C);
	std::cout << "the normal of the plane [" << N <<"]"<< std::endl;
	std::cout << "one point on the plane [" << C << "]" << std::endl;
}

void project_hole_vertex_to_plane(Eigen::MatrixXd &Hole_vertex, Eigen::RowVector3d &N, Eigen::RowVector3d &C,Eigen::MatrixXd &ProjectTo_vertex) {

	ProjectTo_vertex.resize(Hole_vertex.rows(), 3);
	double t; // t = (a*x+b*y+c*z+d)/(a^2+b^2+c^2)
	double a = N(0, 0);
	double b = N(0, 1);
	double c = N(0, 2);
	double x0 = C(0, 0);
	double y0 = C(0, 1);
	double z0 = C(0, 2);
	double d = -1 * (a*x0 + b * y0 + c * z0);
	//double value;
	for (int i = 0; i < Hole_vertex.rows(); i++) {
		t = (a*Hole_vertex(i, 0) + b * Hole_vertex(i, 1) + c * Hole_vertex(i, 2) + d) / (a*a + b * b + c * c);
		ProjectTo_vertex(i, 0) = Hole_vertex(i, 0) -a*t; // x' = x - a * t
		ProjectTo_vertex(i, 1) = Hole_vertex(i, 1) -b*t; // y' = y - b * t
		ProjectTo_vertex(i, 2) = Hole_vertex(i, 2) -c*t; // z' = z - c * t
		
		//// test on the plane? Actually all p are on projected correctly, but due to calculation error, there may be some deviation
		//value = a * ProjectTo_vertex(i, 0) + b * ProjectTo_vertex(i, 1) + c * ProjectTo_vertex(i, 2) + d;

		//if (value!= 0) {
		//	std::cout << "not on the plane" << std::endl;
		//	std::cout << "missing value = " << value << std::endl;
		//}
	}
}

void rotate_to_xy_plane(Eigen::RowVector3d &N, Eigen::MatrixXd &ProjectTo_vertex, Eigen::MatrixXd &vertex_on_xy) {
	Eigen::Vector3d n_plane, n_xy; // n_plane is the normal of the fitted plane, n_xy is the normal of the xy_plane
	Eigen::MatrixXd Rotation_matrix(3, 3); // the rotation matrix
	n_plane = N;
	n_xy = Eigen::Vector3d(0, 0, 1);
	
	// initialize the vertex_on_xy same as Projected_Vertex
	vertex_on_xy.resize(ProjectTo_vertex.rows(),3);
	// calculate the rotation matrix
	Rotation_matrix = igl::rotation_matrix_from_directions(n_plane, n_xy);

	std::cout << "project rows cols " << ProjectTo_vertex.rows() << " " << ProjectTo_vertex.cols() << std::endl;
	std::cout << "v on xy rows cols " << vertex_on_xy.rows() << " " << vertex_on_xy.cols() << std::endl;

	// calculate the rotated vertex coordinate based on the rotation matrix
	// [NOTE] The matrix multiplication requires the transpose operation, 
	// so the vertex_on_xy.transpose() is the final result which same format as normal V from a mesh
	vertex_on_xy = Rotation_matrix * ProjectTo_vertex.transpose();
	std::cout << "vertex on xy " << vertex_on_xy.transpose() << std::endl; 
}

void constrained_delauney_triangulation(Eigen::MatrixXd &vertex_on_xy) {

}

void project_hole_vertex_back() {

}
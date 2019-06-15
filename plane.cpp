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
	double value;
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
	
	//Eigen::Vector3d N_xy = Eigen::RowVector3d(0, 0, 1); // normal of xy plane (0, 0, 1)
	//// test rotation matrix
	//Eigen::Vector3d N_xyy = Eigen::RowVector3d(1, 0, 0); // normal of xy plane (0, 0, 1)
	//Eigen::MatrixXd rotM(3, 3);
	//std::cout << "rotm rows cols " << rotM.rows() << " " << rotM.cols() << std::endl;
	//// initialize the vertex_on_xy matrix
	////vertex_on_xy.resize(ProjectTo_vertex.rows(), 2);
	////std::cout << "vertex-on-xy rows cols " << vertex_on_xy.rows() << " " << vertex_on_xy.cols() << std::endl;
	//igl::rotation_matrix_from_directions(N_xy, N_xyy);

}

void project_hole_vertex_back() {

	Eigen::Vector3d N_xy = Eigen::Vector3d(0, 0, 1); // normal of xy plane (0, 0, 1)
	// test rotation matrix
	Eigen::Vector3d N_xyy = Eigen::Vector3d(1, 0, 0); // normal of xy plane (0, 0, 1)
	std::cout << "vector rows cols " << N_xy.rows() << " " << N_xy.cols() << std::endl;
	Eigen::MatrixXd rotM(3, 3);
	std::cout << "rotm rows cols " << rotM.rows() << " " << rotM.cols() << std::endl;
	// initialize the vertex_on_xy matrix
	//vertex_on_xy.resize(ProjectTo_vertex.rows(), 2);
	//std::cout << "vertex-on-xy rows cols " << vertex_on_xy.rows() << " " << vertex_on_xy.cols() << std::endl;
	rotM = igl::rotation_matrix_from_directions(N_xy, N_xyy);
	std::cout << "rotm " << rotM << std::endl;
}
#include "Header/plane.h"
// test, to be deleted
#include "Header/io.h"

Eigen::MatrixXd vertex_xy_coordinates;
Eigen::MatrixXd t_vertex_on_xy;
//Eigen::MatrixXd bc;

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
	// [NOTE] due to the calc error, the z value may not be exactlly 0 but z will be a constant 
	// which means the rotated vertices are on a plane parallel with xy plane
}


void constrained_delauney_triangulation(Eigen::MatrixXd &vertex_on_xy, Eigen::MatrixXi &cdt_f, Eigen::MatrixXd &bc) {

	// transpose teh vertex_on_xy to get (x, y, z) format and resize based on the rows and 2
	t_vertex_on_xy = vertex_on_xy.transpose();
	vertex_xy_coordinates.resize(t_vertex_on_xy.rows(), 2);

	// slice only col(0)---x and col(1)-----y of vertex_on_xy
	vertex_xy_coordinates.col(0) = t_vertex_on_xy.col(0);
	vertex_xy_coordinates.col(1) = t_vertex_on_xy.col(1);

	std::cout << "v xy coordinate " << vertex_xy_coordinates << std::endl;

	// perform Constrained Delaunay Triangulation
	// define orient2D and incircle functor for DT
	const auto &orient2d_functor = [](const double *pa, const double *pb, const double *pc) {
		return orient2D(pa, pb, pc);
	};
	const auto &incircle_functor = [](const double *pa, const double *pb, const double *pc, const double *pd) {
		return incircle(pa, pb, pc, pd);
	};

	// get basic DT
	igl::delaunay_triangulation(vertex_xy_coordinates, orient2d_functor, incircle_functor, cdt_f); 
	// [NOTE!!!] There seems a bug in igl::lexicographic_triangulation.cpp line 95 to line 108, I've commit on it

	// there are some invalid faces in basic DT which is "out of" the convex hull of the polygon formed by xy coordinates

	// [NOTE] igl::barycenter requires the V be a 3D matrix (x,y,z), add col(2)---z back
	vertex_xy_coordinates.conservativeResize(t_vertex_on_xy.rows(), 3);
	vertex_xy_coordinates.col(2) = t_vertex_on_xy.col(2);
	
	// extract valid delaunay result 
	extract_valid_cdt_f(cdt_f, bc, vertex_xy_coordinates,vertex_xy_coordinates);

	// refinement for the basic delaunay result
	// calculate the average edge length of the poly boundary
	double avg_edge_len;
	for (int i = 0; i < vertex_xy_coordinates.rows(); i++) {
		// TODO
	}
}

bool is_point_in_poly(Eigen::MatrixXd &poly, double &x_bc,double &y_bc) {
	// cast a ray from the query point, count how many edges it crosses
	// if cross is odd-----inside
	// if cross is even----outside
	int crossing = 0;
	double slope;
	bool cond1, cond2, above;

	for (int i = 0; i < poly.rows(); i++) {
		if (i == poly.rows() - 1) {
			// if poly(i) is the last vertex, connect it with the first vertex (x0,y0)
			slope = (poly(0, 1) - poly(i, 1)) / (poly(0, 0) - poly(i, 0));
			cond1 = (poly(i, 0) < x_bc) && (poly(0, 0) > x_bc);
			cond2 = (poly(0, 0) < x_bc) && (poly(i, 0) > x_bc);
		}
		else {
			slope = (poly(i + 1, 1) - poly(i, 1)) / (poly(i + 1, 0) - poly(i, 0));
			cond1 = (poly(i, 0) < x_bc) && (poly(i + 1, 0) > x_bc);
			cond2 = (poly(i + 1, 0) < x_bc) && (poly(i, 0) > x_bc);
		}

		above = (y_bc < slope*(x_bc - poly(i, 0)) + poly(i, 1));
		if ((cond1 || cond2) && above) {
			crossing += 1;
		}
	}
	return (crossing % 2 != 0);
}

void extract_valid_cdt_f(Eigen::MatrixXi &cdt_f, Eigen::MatrixXd &bc, Eigen::MatrixXd &vertex_convex_hull,Eigen::MatrixXd &vertex_all) {
	
	//calculate the barycenter of each face
	bc.resize(cdt_f.rows(), 3);

	igl::barycenter(vertex_all, cdt_f, bc);

	// check if the bc is in the polygon(2D) formed by V, if the bc is outside it means the corresponding face should be deleted
	Eigen::MatrixXd poly_v_2d;
	Eigen::RowVectorXi pos_valid_F; // mark which F is valid
	Eigen::MatrixXi valid_cdt_F; // remove those faces which barycenter is out of the poly area
	int valid_f_count = 0; // how many valid faces
	poly_v_2d.resize(vertex_convex_hull.rows(), 2); // 2D poly
	poly_v_2d.col(0) = vertex_convex_hull.col(0);
	poly_v_2d.col(1) = vertex_convex_hull.col(1);

	color_bc.resize(bc.rows(), 3);
	color_bc.setConstant(1);

	pos_valid_F.resize(cdt_f.rows());
	pos_valid_F.setZero(); // initialize as 0, if 1 means valid face

	for (int i = 0; i < bc.rows(); i++) {
		if (is_point_in_poly(poly_v_2d, bc(i, 0), bc(i, 1))) {
			// if the bc is inside the Poly, keep the face
			std::cout << "#bc " << i << " is inside" << std::endl;
			valid_f_count += 1;
			pos_valid_F(0, i) = 1;
			// visualize the valid barycenter
			color_bc.row(i) << 1, 0, 0;
		}
	}
	std::cout << "valid pos " << pos_valid_F << std::endl;
	// resize according to the valid face count
	valid_cdt_F.resize(valid_f_count, 3);
	std::cout << "#valid faces " << valid_f_count << std::endl;

	int j = 0; // mark the pos in valid_cdt_f
	// extract valid face
	for (int i = 0; i < cdt_f.rows(); i++) {
		if (pos_valid_F(0, i) == 1) {
			std::cout << "keep #f " << i << std::endl;
			valid_cdt_F.row(j) = cdt_f.row(i);
			j += 1;
		}
	}
	cdt_f.resize(valid_f_count, 3);
	cdt_f = valid_cdt_F;
}

void refinement_on_basic_delaunay(Eigen::MatrixXi &cdt_f) {

}
void project_hole_vertex_back() {

}
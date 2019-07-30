#include "Header/plane.h"
// test, to be deleted
#include "Header/io.h"

#define ACCURACY 0.000001
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

void constrained_delauney_triangulation(Eigen::MatrixXd &vertex_on_xy, Eigen::MatrixXi &cdt_f, Eigen::MatrixXd &bc, Eigen::MatrixXd &cdt_v, Eigen::MatrixXd &vertex_new) {

	// transpose teh vertex_on_xy to get (x, y, z) format and resize based on the rows and 2
	t_vertex_on_xy = vertex_on_xy.transpose();
	vertex_xy_coordinates.resize(t_vertex_on_xy.rows(), 2);

	// slice only col(0)---x and col(1)-----y of vertex_on_xy
	vertex_xy_coordinates.col(0) = t_vertex_on_xy.col(0);
	vertex_xy_coordinates.col(1) = t_vertex_on_xy.col(1);

	//std::cout << "v xy coordinate " << vertex_xy_coordinates << std::endl;

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
	
	
	Eigen::MatrixXd dblA; // all subtriangle's [doublearea], the subtriangle's area is 1/2 of each element

	// calculate each triangle's area, send dblA is to avoid thin triangle
	igl::doublearea(vertex_xy_coordinates, cdt_f, dblA); 

	// extract valid delaunay result 
	extract_valid_cdt_f(cdt_f, bc, vertex_xy_coordinates,vertex_xy_coordinates,dblA);

	// update the dblA for extracted valid faces
	igl::doublearea(vertex_xy_coordinates, cdt_f, dblA);

	// refinement for the basic delaunay result
	// calculate the average edge length of the poly boundary
	double convex_length = 0;
	for (int i = 0; i < vertex_xy_coordinates.rows(); i++) {
		if (i == vertex_xy_coordinates.rows() - 1) {
			// the last edge is v_last & v0
			convex_length += (vertex_xy_coordinates.row(i) - vertex_xy_coordinates.row(0)).norm();
		}
		else {
			convex_length += (vertex_xy_coordinates.row(i) - vertex_xy_coordinates.row(i + 1)).norm();
		}
	}
	convex_length = convex_length / (double)vertex_xy_coordinates.rows();
	std::cout << "avg convex len = " << convex_length << std::endl;
	double epsilon = (sqrt(3)/4)*convex_length*convex_length; // ¦Å=(¡Ì3/4)*avg_len^2
	std::cout << "epsilon = " << epsilon << std::endl;

	// perform Constrained Delaunay Refinement, the constraint is that any subtriangle's area <= ¦Å
	Eigen::MatrixXd v_of_subtriangle(3, 3); // the three v of the triangle
	//Eigen::MatrixXd vertex_new; // store newly added vertex, it will finally append to the vertex_xy_coordinates
	int num_new_v = 0; // count how many new v are created
	for (int i = 0; i < cdt_f.rows(); i++) {
		if (dblA(i, 0) > epsilon*2) {
			// the subtriangle is large than ¦Å, needs to be subdivided
			//std::cout << "#triangle " << i << " is large. " << std::endl;
			v_of_subtriangle.row(0) = vertex_xy_coordinates.row(cdt_f(i, 0));
			v_of_subtriangle.row(1) = vertex_xy_coordinates.row(cdt_f(i, 1));
			v_of_subtriangle.row(2) = vertex_xy_coordinates.row(cdt_f(i, 2));
			refinement_on_basic_delaunay(v_of_subtriangle, epsilon, vertex_new, num_new_v);
		}
	}

	cdt_v.resize(vertex_xy_coordinates.rows() + vertex_new.rows(), 2);
	cdt_v <<
		vertex_xy_coordinates.col(0), vertex_xy_coordinates.col(1),
		vertex_new.col(0), vertex_new.col(1);

	// perform another delaunay triangulation
	igl::delaunay_triangulation(cdt_v, orient2d_functor, incircle_functor, cdt_f);
	cdt_v.conservativeResize(cdt_v.rows(), 3);
	cdt_v.col(2).setConstant(vertex_new(0, 2));

	// extract valid faces
	igl::doublearea(cdt_v, cdt_f, dblA);
	extract_valid_cdt_f(cdt_f, bc, vertex_xy_coordinates, cdt_v, dblA);
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

void extract_valid_cdt_f(Eigen::MatrixXi &cdt_f, Eigen::MatrixXd &bc, Eigen::MatrixXd &vertex_convex_hull,Eigen::MatrixXd &vertex_all, Eigen::MatrixXd &dblA) {
	
	//calculate the barycenter of each face
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
		if (is_point_in_poly(poly_v_2d, bc(i, 0), bc(i, 1))&&((dblA(i, 0)*0.5) > ACCURACY)) {
			// if the bc is inside the Poly and not a 'thin triangle', keep the face
			//std::cout << "#bc " << i << " is inside" << std::endl;
			valid_f_count += 1;
			pos_valid_F(0, i) = 1;
			// visualize the valid barycenter
			color_bc.row(i) << 1, 0, 0;
		}
	}
	//std::cout << "valid pos " << pos_valid_F << std::endl;
	// resize according to the valid face count
	valid_cdt_F.resize(valid_f_count, 3);
	//std::cout << "#valid faces " << valid_f_count << std::endl;

	int j = 0; // mark the pos in valid_cdt_f
	// extract valid face
	for (int i = 0; i < cdt_f.rows(); i++) {
		if (pos_valid_F(0, i) == 1) {
			//std::cout << "keep #f " << i << std::endl;
			valid_cdt_F.row(j) = cdt_f.row(i);
			j += 1;
		}
	}
	cdt_f.resize(valid_f_count, 3);
	cdt_f = valid_cdt_F;
}

void refinement_on_basic_delaunay(Eigen::MatrixXd &vertex_of_triangle, double &epsilon, Eigen::MatrixXd &vertex_append, int &num_v) {
	// add a new vertex at the barycenter of the triangle
	num_v += 1;
	Eigen::MatrixXd sub_bc;
	Eigen::MatrixXd check_new_subtriangle_v(4, 3);
	Eigen::MatrixXi check_new_subtriangle_f(3, 3);
	igl::barycenter(vertex_of_triangle, Eigen::RowVector3i(0, 1, 2), sub_bc);
	// add the barycenter as a new point
	vertex_append.conservativeResize(num_v, 3);
	vertex_append.row(num_v - 1) = sub_bc;

	// check the new created 3 triangles
	check_new_subtriangle_v.row(0) = vertex_of_triangle.row(0);
	check_new_subtriangle_v.row(1) = vertex_of_triangle.row(1);
	check_new_subtriangle_v.row(2) = vertex_of_triangle.row(2);
	check_new_subtriangle_v.row(3) = sub_bc;
	check_new_subtriangle_f <<
		3, 0, 1,
		3, 1, 2,
		3, 2, 0;

	//std::cout << "vertex of tri " << vertex_of_triangle << std::endl;
	//std::cout << "bc " << sub_bc << std::endl;
	//std::cout << "append new vertex [ " << vertex_append.row(num_v - 1) << " ]" << std::endl;
	//std::cout << "check new subtri v " << check_new_subtriangle_v << std::endl;
	//std::cout << "current append vertex matrix " << vertex_append << std::endl;

	// check the area of the three 
	Eigen::MatrixXd sub_dblA; // store the area of each sub triangle
	Eigen::MatrixXd temp_v(3, 3); // make a proper parameter for refinment_on_basic_delaunay()
	igl::doublearea(check_new_subtriangle_v, check_new_subtriangle_f, sub_dblA);
	//std::cout << "sub dblA " << sub_dblA << std::endl;
	for (int i = 0; i < 3; i++) {
		if (sub_dblA(i, 0) > epsilon * 2) {
			//std::cout << "need further process" << std::endl;
			temp_v.row(0) = check_new_subtriangle_v.row(check_new_subtriangle_f(i, 0));
			temp_v.row(1) = check_new_subtriangle_v.row(check_new_subtriangle_f(i, 1));
			temp_v.row(2) = check_new_subtriangle_v.row(check_new_subtriangle_f(i, 2));
			refinement_on_basic_delaunay(temp_v, epsilon, vertex_append, num_v);
		}
	}
}

void project_hole_vertex_back(Eigen::MatrixXd &cdt_vertex, Eigen::MatrixXi &cdt_face, Eigen::MatrixXd &v_original_boundary, Eigen::MatrixXd &vertex_new, Eigen::MatrixXd &vertex_new_3D) {
	// calculate the adjacent list of each vertex
	std::vector< std::vector<double> > adj; // store the adjacent info, adj[i][] means the adjacent vertex of v(i)
	igl::adjacency_list(cdt_face,adj,true); 

	// Mean value Coordinates
	Eigen::RowVectorXd wi, ni, tan_half_value, internal_edge_len; // wi/ni: for each planar distribution
	Eigen::MatrixXd ni_total; // reconstructed matrix of all ni for computing
	ni_total.resize(vertex_new.rows(), cdt_vertex.rows());
	ni_total.setZero();


	int count_ni_total = 0; // control which row of ni_total to construct
	int count_ni; // control which col of ni to construct


	double a, b, c, sum_wi;
	int num_v_orignial_boundary = v_original_boundary.rows();

	//for (int i = num_v_orignial_boundary; i < adj.size(); i++) {
	//	for (int j = 0; j < adj[i].size(); j++) {
	//		std::cout << adj[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}


	for (int i = num_v_orignial_boundary; i < adj.size(); i++) {
		// v before num_v_o_b is on original mesh boundary, use there 3D coordinates directly 
		sum_wi = 0;
		wi.resize(adj[i].size());
		ni.resize(adj[i].size());
		tan_half_value.resize(wi.cols());
		internal_edge_len.resize(wi.cols());


		for (int j = 0; j < adj[i].size(); j++) {

			// for each adjacent v, tan a/2 is v_o---v0---v1
			if (j == adj[i].size() - 1) {
				// for the last adjacent v, edge a is vj---v0
				a = sqrt((cdt_vertex(adj[i][j], 0) - cdt_vertex(adj[i][0], 0))*(cdt_vertex(adj[i][j], 0) - cdt_vertex(adj[i][0], 0)) +
					(cdt_vertex(adj[i][j], 1) - cdt_vertex(adj[i][0], 1))*(cdt_vertex(adj[i][j], 1) - cdt_vertex(adj[i][0], 1)));
			}
			else {
				a = sqrt((cdt_vertex(adj[i][j], 0) - cdt_vertex(adj[i][j + 1], 0))*(cdt_vertex(adj[i][j], 0) - cdt_vertex(adj[i][j + 1], 0)) +
					(cdt_vertex(adj[i][j], 1) - cdt_vertex(adj[i][j + 1], 1))*(cdt_vertex(adj[i][j], 1) - cdt_vertex(adj[i][j + 1], 1)));

			}

			if (j == 0) {
				// for the first adjacent v, the edge b needs to be calc
				b = sqrt((cdt_vertex(adj[i][j], 0) - cdt_vertex(i, 0))*(cdt_vertex(adj[i][j], 0) - cdt_vertex(i, 0)) +
					(cdt_vertex(adj[i][j], 1) - cdt_vertex(i, 1))*(cdt_vertex(adj[i][j], 1) - cdt_vertex(i, 1)));

				internal_edge_len(0, j) = b;

			}
			else {
				// for other adjacent v, edge b has already been calculated in last loop
				b = internal_edge_len(0, j);
			}

			if (j == adj[i].size() - 1) {
				// for the last adjacent v, the edge c has been calculated in the first loop
				c = internal_edge_len(0, 0);
			}
			else {
				c = sqrt((cdt_vertex(i, 0) - cdt_vertex(adj[i][j + 1], 0))*(cdt_vertex(i, 0) - cdt_vertex(adj[i][j + 1], 0)) +
					(cdt_vertex(i, 1) - cdt_vertex(adj[i][j + 1], 1))*(cdt_vertex(i, 1) - cdt_vertex(adj[i][j + 1], 1)));

				internal_edge_len(0, j + 1) = c;
			}
			// calc tan a/2
			tan_half_value(0, j) = tan_half_angle(a, b, c);

			// as for tan(a-1)/2, only need to calc for the first adjacent v
			if(j==0){
				// for the first adjacent v, tan (a-1)/2 is v_o---v_max---v0
				a = sqrt((cdt_vertex(adj[i][0], 0) - cdt_vertex(adj[i][adj[i].size() - 1], 0))*(cdt_vertex(adj[i][0], 0) - cdt_vertex(adj[i][adj[i].size() - 1], 0)) +
					(cdt_vertex(adj[i][0], 1) - cdt_vertex(adj[i][adj[i].size() - 1], 1))*(cdt_vertex(adj[i][0], 1) - cdt_vertex(adj[i][adj[i].size() - 1], 1)));
				b = sqrt((cdt_vertex(adj[i][adj[i].size() - 1], 0) - cdt_vertex(i, 0))*(cdt_vertex(adj[i][adj[i].size() - 1], 0) - cdt_vertex(i, 0)) +
					(cdt_vertex(adj[i][adj[i].size() - 1], 1) - cdt_vertex(i, 1))*(cdt_vertex(adj[i][adj[i].size() - 1], 1) - cdt_vertex(i, 1)));
				c = internal_edge_len(0, j);

				internal_edge_len(0, adj[i].size() - 1) = b;
				tan_half_value(0, adj[i].size() - 1) = tan_half_angle(a, b, c);
				// calc wi, wi = (tan(a-1)/2+tana/2)/(||vi-vo||)
				wi(0, j) = (tan_half_value(0, adj[i].size() - 1) + tan_half_value(0, j)) / internal_edge_len(0, j);
				sum_wi += wi(0, j);
			}
			else {
				// for other adjacent v, tan(a-1)/2 has already been calculated in previous loop when calc tan a/2
				// calc wi, wi = (tan(a-1)/2+tana/2)/(||vi-vo||)
				wi(0, j) = (tan_half_value(0, j - 1) + tan_half_value(0, j)) / internal_edge_len(0, j);
				sum_wi += wi(0, j);
			}
		}

		//std::cout << "current wi: " << wi << std::endl;

		// [TODO] calculate ni
		for (int i = 0; i < wi.cols(); i++) {
			
			ni(0, i) = wi(0, i) / sum_wi;
		}
		//std::cout << "current ni: " << ni << std::endl;
		
		// construct the ni_total
		for (int k = 0; k < ni.cols(); k++) {
			// update actual ni value for marked position
			ni_total(count_ni_total, adj[i][k]) = ni(0, k);
		}

		count_ni_total += 1;
	}

	// [TODO]solve the equaltion to get 3D vertices coordinates
	// construct Ax = b
	Eigen::MatrixXd A,B, B_left,A_left, A_right;

	B_left = ni_total.leftCols(num_v_orignial_boundary);
	B = B_left * v_original_boundary;

	A_left.resize(count_ni_total, count_ni_total);
	A_left.setIdentity();
	A_right = ni_total.rightCols(count_ni_total);
	A = A_left - A_right;


	vertex_new_3D = A.colPivHouseholderQr().solve(B);

}

double tan_half_angle(double &len_a, double &len_b, double &len_c) {
	double tan_half;
	tan_half = sqrt(((len_a - len_b + len_c)*(len_a + len_b - len_c)) / ((len_a + len_b + len_c)*(-len_a + len_b + len_c)));
	return tan_half;
}

void seampatch(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &v_on_line_R, Eigen::MatrixXd &v_new_3D_R,
	Eigen::MatrixXi &cdt_face_R, Eigen::RowVectorXi &hole_idx_R, Eigen::MatrixXd &v_on_line_L, Eigen::MatrixXd &v_new_3D_L,
	Eigen::MatrixXi &cdt_face_L, Eigen::RowVectorXi &hole_idx_L) {

	Eigen::MatrixXd V_total;
	Eigen::MatrixXi F_total;

	int numv_ol_R = v_on_line_R.rows();
	int numv_ol_L = v_on_line_L.rows();
	int numv_R = v_new_3D_R.rows();
	int numv_L = v_new_3D_L.rows();
	int numv_O = V.rows();
	int numf_R = cdt_face_R.rows();
	int numf_L = cdt_face_L.rows();
	int numf_O = F.rows();

	std::cout << "count covert left = " << count_cover_hole_part_left << std::endl;
	std::cout << "count covert right = " << count_cover_hole_part_right << std::endl;

	std::cout << "count new = " << count_cover_new << std::endl;
	// combine all vertices, V + v_on_line_right + v_on_line_left + new_v_right + new_v_left
	// combine all faces, F + new_f_right + new_f_left
	V_total.resize(numv_O + numv_ol_R + numv_R + numv_ol_L + numv_L, 3);
	F_total.resize(numf_O + numf_R + numf_L, 3);

	V_total <<
		V,
		v_on_line_R,
		v_new_3D_R,
		v_on_line_L,
		v_new_3D_L;

	// adjust the idx of cdt faces(vertex's idx)
	for (int i = 0; i < numf_R; i++) {
		for (int j = 0; j < 3;j++) {
			if (cdt_face_R(i, j) < hole_idx_R.cols()) {
				// if #cdt_f idx < # hole, it means the cdt_f vertex is boundary vertex
				cdt_face_R(i, j) = hole_idx_R(0, cdt_face_R(i, j));
			}
			else {
				// if #cdt_f idx >= # hole, it means the cdt_f vertex is either v_online or new_v
				cdt_face_R(i, j) = cdt_face_R(i, j) - hole_idx_R.cols() + numv_O;
			}
		}
	}

	std::cout << "cdt face l " << cdt_face_L << std::endl;
	switch (cover_origin)
	{
	case 0:
		// if the idx of selected vertex not cover the origin point
		for (int i = 0; i < numf_L; i++) {
			for (int j = 0; j < cdt_face_L.cols(); j++) {
				if (cdt_face_L(i, j) < hole_idx_L.cols()) {
					// if #cdt_f idx < # hole, it means the cdt_f vertex is boundary vertex
					cdt_face_L(i, j) = hole_idx_L(0, cdt_face_L(i, j));
				}
				else {
					cdt_face_L(i, j) = cdt_face_L(i, j) - hole_idx_L.cols() + numv_O + numv_ol_R + numv_R;
				}
			}
		}
		break;
	case 1:
		// if the idx of selected vertex not cover the origin point
		for (int i = 0; i < numf_L; i++) {
			for (int j = 0; j < cdt_face_L.cols(); j++) {
				if (cdt_face_L(i, j) < count_cover_hole_part_left) {
					// if #cdt_f idx < # hole, it means the cdt_f vertex is boundary vertex
					cdt_face_L(i, j) = hole_idx_L(0, cdt_face_L(i, j));
				}
				else if (cdt_face_L(i, j) >= count_cover_hole_part_left && cdt_face_L(i, j) < (count_cover_hole_part_left + count_cover_new)) {
					cdt_face_L(i, j) = cdt_face_L(i, j) - count_cover_hole_part_left + numv_O + numv_ol_R + numv_R;

				}
				else if (cdt_face_L(i, j) >= (count_cover_hole_part_left + count_cover_new) && cdt_face_L(i, j) < (count_cover_hole_part_left + count_cover_new + count_cover_hole_part_right)) {
					cdt_face_L(i, j) = hole_idx_L(cdt_face_L(i, j) - count_cover_new);
				}
				else {
					cdt_face_L(i, j) = cdt_face_L(i, j) - count_cover_hole_part_left - count_cover_hole_part_right + numv_O + numv_ol_R + numv_R;
				}
			}
		}
		break;
	default:
		break;
	}



	F_total <<
		F,
		cdt_face_R,
		cdt_face_L;


	// make orientation of cdt face same as original mesh
	Eigen::MatrixXi F_orient;
	igl::bfs_orient(F_total, F_orient, V);


	V.resize(numv_O + numv_ol_R + numv_R + numv_ol_L + numv_L, 3);
	F.resize(numf_O + numf_R + numf_L, 3);
	V = V_total;
	//F = F_total;
	F = F_orient;


}


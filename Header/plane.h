#pragma once
/*
	Created by Yuhan Ping from HKU on 2019-06-13.
	This is the header file for fitting plane and project the hole vertex from 3D to 2D.
	Under development.
	Latest updated on 2019-07-08.
*/
#ifndef PLANE_H
#define PLANE_H

#include <vector>
#include <igl/fit_plane.h>
#include <igl/rotation_matrix_from_directions.h>
#include <igl/delaunay_triangulation.h>
#include <igl/barycenter.h>
#include <igl/point_in_poly.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/edge_lengths.h>
#include <igl/doublearea.h>
#include <igl/adjacency_list.h>


#include <iostream>
#include <fstream>





void get_plane(Eigen::MatrixXd &Hole_vertex, Eigen::RowVector3d &N, Eigen::RowVector3d &C);
/*
	fit a plane for hole vertices
	[Return value] 1) a normal of the plane; 2) a point on the plane
*/

void project_hole_vertex_to_plane(Eigen::MatrixXd &Hole_vertex, Eigen::RowVector3d &N, Eigen::RowVector3d &C, Eigen::MatrixXd &ProjectTo_vertex);
/*
	project the hole vertex to the plane
	[Return value] Projected_vertex
*/

void rotate_to_xy_plane(Eigen::RowVector3d &N, Eigen::MatrixXd &ProjectTo_vertex, Eigen::MatrixXd &vertex_on_xy);
/*
	rotate the plane to a xy plane for further Constrained Delauney Triangulation
	[Return value] vertex on the xy plane
*/

void constrained_delauney_triangulation(Eigen::MatrixXd &vertex_on_xy, Eigen::MatrixXi &cdt_f,Eigen::MatrixXd &bc, Eigen::MatrixXd &cdt_v, Eigen::MatrixXd &vertex_new);
/*
	Constrained Delauney Triangulation on the xy plane 
	[Return value] triangulated vertex and face
*/

template <typename T>
inline int orient2D(const T &pa, const T &pb, const T &pc) {
	/*
	A functor such that orient2D(pa, pb, pc) returns
				[]1 if pa,pb,pc forms a conterclockwise triangle.
				[-1] if pa,pb,pc forms a clockwise triangle.
				[0] if pa,pb,pc are collinear.
				where the argument pa,pb,pc are of type Scalar[2].
	*/
	// pa & pb & pc : (x, y)
	// calculate the slope of v(ab) and v(ac)
	double vab_x = pb[0] - pa[0];
	double vab_y = pb[1] - pa[1];
	double vbc_x = pc[0] - pb[0];
	double vbc_y = pc[1] - pb[1];
	double check_slope = vab_y * vbc_x - vbc_y * vab_x;
	if (check_slope < 0) {
		// slope vab < slope vbc, counterclockwise
		//std::cout << "counter clockwise" << std::endl;
		return 1;
	}
	else if (check_slope > 0) {
		// slope vab > slope vbc, clockwise
		//std::cout << "clockwise" << std::endl;

		return -1;
	}
	else {
		// slope vab = slope vbc, collinear
		//std::cout << "collinear" << std::endl;

		return 0;
	}
	return 0;
}


template <typename T>
int incircle(const T &pa, const T &pb, const T &pc,const T &pd) {
	/*
	A functor such that incircle(pa, pb, pc, pd) returns
				 [1] if pd is on the positive size of circumcirle of (pa,pb,pc)
				 [-1] if pd is on the positive size of circumcirle of (pa,pb,pc)
				 [0] if pd is cocircular with pa, pb, pc.
	*/
	double vad_x = pa[0] - pd[0];
	double vad_y = pa[1] - pd[1];
	double vbd_x = pb[0] - pd[0];
	double vbd_y = pb[1] - pd[1];
	double vcd_x = pc[0] - pd[0];
	double vcd_y = pc[1] - pd[1];
	double angle_abd = vad_x * vbd_y - vbd_x * vad_y;
	double angle_bcd = vbd_x * vcd_y - vcd_x * vbd_y;
	double angle_cad = vcd_x * vad_y - vad_x * vcd_y;
	double len_ad = vad_x * vad_x + vad_y * vad_y;
	double len_bd = vbd_x * vbd_x + vbd_y * vbd_y;
	double len_cd = vcd_x * vcd_x + vcd_y * vcd_y;
	double check_angle = len_ad * angle_bcd + len_bd * angle_cad + len_cd * angle_abd;
	if (check_angle < 0) {
		// outside the circumcircle
		//std::cout << "outside" << std::endl;
		return -1;
	}
	else if (check_angle > 0) {
		// inside the circle
		//std::cout << "inside" << std::endl;
		return 1;
	}
	else {
		// on the circle
		//std::cout << "cocircular" << std::endl;
		return 0;
	}
	return 0;
}

bool is_point_in_poly(Eigen::MatrixXd &poly,double &x_bc, double &y_bc);
/*
	check if the barycenter is inside the Polygon by V (2D)
	[Return Value] True---inside, false---outside
*/

void extract_valid_cdt_f(Eigen::MatrixXi &cdt_f, Eigen::MatrixXd &bc, Eigen::MatrixXd &vertex_convex_hull, Eigen::MatrixXd &vertex_all, Eigen::MatrixXd &dblA);
/*
	extract valid delaunay faces
	[Return value] valid delaunay faces
*/

void refinement_on_basic_delaunay(Eigen::MatrixXd &vertex_of_triangle, double &epsilon, Eigen::MatrixXd &vertex_append, int &num_v);
/*
	perform refinment on basic Delaunay result following the constraint that any subtriangle's area <= ¦Å
	¦Å=(¡Ì3/4)*avg_len^2
	[Return Value] the final delaunay faces which follow specific constraints
*/

void project_hole_vertex_back(Eigen::MatrixXd &cdt_vertex, Eigen::MatrixXi &cdt_face, Eigen::MatrixXd &v_original_boundary, Eigen::MatrixXd &vertex_new_2D, Eigen::MatrixXd &vertex_new_3D);
/*
	use 'mean value coordinates' to project the CDT vertices back to 3D coordinates
	[Return value] 3D coordinates (x,y,z) of the CDT result vertices
*/

double tan_half_angle(double &len_a, double &len_b, double &len_c);
/*
	calculate the tangent of half of angle BAC, len_a/b/c is the length of the edge of the triangle
	[Return value] tan a/2
*/


#endif // !PLANE_H


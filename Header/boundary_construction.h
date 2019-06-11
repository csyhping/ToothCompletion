#pragma once
/*
	Created by Yuhan Ping from HKU on 2019-06-06.
	This is the header file for boundary construction.
	Under development.
	Latest updated on 2019-06-11.
*/
#ifndef BOUNDARY_CONSTRUCTION_H
#define BOUNDARY_CONSTRUCTION_H

#include <igl/boundary_loop.h>
#include <igl/avg_edge_length.h>
#include <igl/edge_lengths.h>
#include <igl/boundary_loop.h>

extern Eigen::MatrixXd Color_per_vertex; // color of each vertex, RED for boundary vertex
extern Eigen::MatrixXd Hole_vertex; // hole_boundary vertices, including select_v1 + select_v2 + new_v_on_line + orginal_boundary_v_above_v1&v2


double calc_average_edge_length(Eigen::MatrixXd &V, Eigen::MatrixXi &F); 
/*
	calculate the average length between two adjacent vertices
	[Return value] avg_length: the average length 
*/

void create_vertex_on_line(Eigen::RowVector3d &select_v1, Eigen::RowVector3d &select_v2, Eigen::MatrixXd &New_v_on_line);
/*
	// create new vertices on the interactive line with the average length 
	[Return value] new vertices coordinates stored in a matrix
*/

void get_hole_boundary(Eigen::MatrixXd &color_per_v, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::RowVector3d &select_v1, Eigen::RowVector3d &select_v2, Eigen::MatrixXd &New_v_on_line);
/*
	// hole_boundary = select_v1 + select_v2 + new_created_v_on_straight_line + boundary_v_above_v1_and_v2
	[Return value] constructed hole boundary vertices
*/

#endif // !BOUNDARY_CONSTRUCTION_H


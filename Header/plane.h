#pragma once
/*
	Created by Yuhan Ping from HKU on 2019-06-13.
	This is the header file for fitting plane and project the hole vertex from 3D to 2D.
	Under development.
	Latest updated on 2019-06-16.
*/
#ifndef PLANE_H
#define PLANE_H
#include <igl/fit_plane.h>
#include <igl/rotation_matrix_from_directions.h>
#include <igl/delaunay_triangulation.h>



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

void constrained_delauney_triangulation();
/*
	Constrained Delauney Triangulation on the xy plane 
	[Return value] triangulated vertex and face
*/

void project_hole_vertex_back();

#endif // !PLANE_H


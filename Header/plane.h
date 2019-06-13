#pragma once
/*
	Created by Yuhan Ping from HKU on 2019-06-13.
	This is the header file for fitting plane and project the hole vertex from 3D to 2D.
	Under development.
	Latest updated on 2019-06-13.
*/
#ifndef PLANE_H
#define PLANE_H
#include <igl/fit_plane.h>


void get_plane(Eigen::MatrixXd &Hole_vertex, Eigen::RowVector3d &N, Eigen::RowVector3d &C);
/*
	fit a plane for hole vertices
	[Return value] 1) a normal of the plane; 2) a point on the plane
*/
void project_hole_vertex_to_plane();
/*
	project the hole vertex to the plane
*/
void project_hole_vertex_back();

#endif // !PLANE_H


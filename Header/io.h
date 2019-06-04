#pragma once
/*
	Created by Yuhan Ping from HKU on 2019-06-04.
	This is the header file for all io-related functions and the viewer settings.
	Under development.
	Latest updated on 2019-06-04.
*/

#ifndef IO_H
#define IO_H

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/unproject_onto_mesh.h>

extern Eigen::Vector3f select_v1, select_v2; // select two vertex to create interactive straight line
extern Eigen::MatrixXd V1;
extern Eigen::MatrixXi F1;


bool mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier); // mouse interaction to select vertex,Eigen::Vector3f &v1, Eigen::Vector3f &v2

bool viewer_display(Eigen::MatrixXd &V, Eigen::MatrixXi &F); // viewer settings

bool load_mesh(std::string &mesh_dir, Eigen::MatrixXd &V, Eigen::MatrixXi &F); // load mesh 


#endif // __IO_H__

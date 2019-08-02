#pragma once
/*
	Created by Yuhan Ping from HKU on 2019-06-04.
	This is the header file for all io-related functions and the viewer settings.
	Under development.
	Latest updated on 2019-07-30.
*/

#ifndef IO_H
#define IO_H

#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/writeOFF.h>
#include <igl/writeDMAT.h>



extern Eigen::MatrixXd V1; 
extern Eigen::MatrixXi F1;
extern int cover_origin;// flag: if the selected two points v3 and v4 cover the origin point of boudnary vertex 
extern int count_cover_hole_part_left, count_cover_new, count_cover_hole_part_right; // mark the position for adjusting cdt_f coordinates when seampatch

// test for visualization, will be deteled finally
extern Eigen::MatrixXd color_bc;
extern std::string inputmesh, prefair_R, prefair_L, postfair_X, prefair_file, postfair_file,hole_R,hole_L;


bool mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier); // mouse interaction to select vertex,Eigen::Vector3f &v1, Eigen::Vector3f &v2

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier); // key interaction to further process after select two vertex

bool viewer_display(Eigen::MatrixXd &V, Eigen::MatrixXi &F); // viewer settings

bool load_mesh(std::string &mesh_dir, Eigen::MatrixXd &V, Eigen::MatrixXi &F); // load mesh 


#endif // __IO_H__

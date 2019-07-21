#pragma once
/*
	Created by Yuhan Ping from HKU on 2019-07-10.
	This is the header file for fairing.
	Under development.
	Latest updated on 2019-07-21.
*/
#ifndef FAIRING_H
#define FAIRING_H

#include <vector>
#include <igl/adjacency_list.h>


#include <iostream>
#include <fstream>



void mesh_fairing(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &one_ring, Eigen::RowVectorXi &hole_idx_R, Eigen::RowVectorXi &hole_idx_L, int &num_new);
/*
	fair the patched area to make it smoother with original mesh
	[Return value] the mesh after fairing
*/

void get_onering(Eigen::RowVectorXi &hole_idx_R, Eigen::RowVectorXi &hole_idx_L, std::vector< std::vector<double> > &adjlist, Eigen::RowVectorXi &posmark, int &num_new, int &num_total);
/*
	get the one-ring vertex of hole vertex
	[Return value] one-ring vertex
*/
#endif // !FAIRING_H


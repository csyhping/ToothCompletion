#include "Header/fairing.h"
#include "Header/io.h"





void mesh_fairing(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &one_ring, Eigen::RowVectorXi &hole_idx_R, Eigen::RowVectorXi &hole_idx_L, int &num_new) {
	// get the adjacent list of each vertex
	std::vector< std::vector<double> > adjlist;
	igl::adjacency_list(F, adjlist, true);

	Eigen::RowVectorXi posmark; // mark which vertex ifx is needed
	posmark.resize(adjlist.size());
	posmark.setZero();

	Eigen::MatrixXd one_LP, two_LP; // laplacian operator. second-order laplacian operator
	Eigen::MatrixXd wi; // weight matrix
	Eigen::MatrixXd len; // length of edge, scale-based wi = ||vi - vj||
	int num_total = 0; // hole + one-ring + new_Vertex

	Eigen::RowVectorXi pos; // the match between original idx and constructed idx
	pos.resize(V.rows());
	pos.setConstant(-1);

	// get one-ring vertex of hole vertex
	get_onering(hole_idx_R, hole_idx_L, adjlist, posmark, num_new, num_total);

	// vertex on line & new vertex, mark as 3
	for (int i = num_new; i < V.rows(); i++) {
		posmark(i) = 3;
		num_total += 1;
	}

	// initialize wi and len
	wi.resize(V.rows() - num_new, num_total);
	wi.setZero();
	len.resize(num_total, num_total); // len(i, j) = len||vi - vj||

	// match between original idx and constructed idx
	int position = 0;
	for (int i = 0; i < posmark.cols(); i++) {
		if (posmark(i) != 0) {
			pos(i) = position;
			position += 1;
		}
	}
	std::cout << "pos " << pos << std::endl;
	std::cout << "v rows cols " << V.rows() << " " << V.cols() << std::endl;
	// calculate w(vi, vj)
	for (int i = 0; i < adjlist.size(); i++) {
		if (posmark(i) == 1) {
			// only focus on 1-ring, boundary and new vertex
			std::cout << "#idx = " << i << "at " << pos(i) << std::endl;

			for (int j = 0; j < adjlist[i].size(); j++) {
				std::cout << "adj[][] = " << adjlist[i][j] << std::endl;
				std::cout << "pos(i) " << pos(i) << std::endl;
				std::cout << "pos(adj) " << pos(adjlist[i][j]) << std::endl;
				std::cout << "v (i,0) " << V(i, 0) << std::endl;
				std::cout << "v(adj) " << V(adjlist[i][j], 0) << std::endl;
				len(pos(i), pos(adjlist[i][j])) = sqrt(
					(V(i, 0) - V(adjlist[i][j], 0))*(V(i, 0) - V(adjlist[i][j], 0)) +
					(V(i, 1) - V(adjlist[i][j], 1))*(V(i, 1) - V(adjlist[i][j], 1)) +
					(V(i, 2) - V(adjlist[i][j], 2))*(V(i, 2) - V(adjlist[i][j], 2))
				);
				len(pos(adjlist[i][j]), pos(i)) = len(pos(i), pos(adjlist[i][j]));
			}
		}
		//getchar();
	}


}


void get_onering(Eigen::RowVectorXi &hole_idx_R, Eigen::RowVectorXi &hole_idx_L, std::vector< std::vector<double> > &adjlist, Eigen::RowVectorXi &posmark, int &num_new, int &num_total) {
	// mark hole idx L, mark as 1
	for (int i = 0; i < hole_idx_L.cols(); i++) {
		posmark(hole_idx_L(i)) = 1;
		num_total += 1;
	}

	// mark hole idx R, mark as 1
	for (int i = 0; i < hole_idx_R.cols(); i++) {
		posmark(hole_idx_R(i)) = 1;
		num_total += 1;
	}

	Eigen::RowVectorXi hole_idx; // collection of left and right hole idx
	hole_idx.resize(hole_idx_L.cols() + hole_idx_R.cols());
	hole_idx <<
		hole_idx_R, hole_idx_L;

	for (int i = 0; i < hole_idx.cols(); i++) {
		for (int j = 0; j < adjlist[hole_idx(i)].size(); j++) {
			if (adjlist[hole_idx(i)][j] < num_new && posmark(adjlist[hole_idx(i)][j]) != 1&& posmark(adjlist[hole_idx(i)][j]) != 2) {
				// #one-ring idx is in original mesh, mark as 2
				posmark(adjlist[hole_idx(i)][j]) = 2;
				num_total += 1;
			}
		}
	}
}
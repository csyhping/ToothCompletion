#include "Header/fairing.h"
#include "Header/io.h"





void mesh_fairing(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::RowVectorXi &hole_idx_R, Eigen::RowVectorXi &hole_idx_L, int &num_new) {
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
	int count_new = 0; // # new vertex
	Eigen::RowVectorXd sum_wi; // sum of wi
	Eigen::MatrixXd A, X, B;

	Eigen::RowVectorXi pos, pos_back; // the match between original idx and constructed idx
	pos.resize(V.rows());
	pos.setConstant(-1);

	// get one-ring vertex of hole vertex
	get_onering(hole_idx_R, hole_idx_L, adjlist, posmark, num_new, num_total);

	// vertex on line & new vertex, mark as 3
	for (int i = num_new; i < V.rows(); i++) {
		posmark(i) = 3;
		num_total += 1;
		count_new += 1;
	}

	// initialize wi and len 
	wi.resize(V.rows() - num_new, num_total);
	wi.setZero();
	len.resize(num_total, num_total); // len(i, j) = len||vi - vj||
	len.setZero();
	sum_wi.resize(num_total);
	sum_wi.setZero();
	pos_back.resize(num_total);

	// initialize AX=B
	A.resize(count_new, count_new);
	X.resize(count_new, 3);
	B.resize(count_new, 3);


	// match between original idx and constructed idx
	int position = 0;
	for (int i = 0; i < posmark.cols(); i++) {
		if (posmark(i) !=0) {
			pos(i) = position;
			pos_back(position) = i;
			position += 1;
		}
	}

	// calculate w(vi, vj)
	for (int i = 0; i < adjlist.size(); i++) {
		if (posmark(i) == 1 || posmark(i) == 3) {
			// only focus on boundary and new vertex, one-ring vertex is only for calc wi of boundary vertex
			for (int j = 0; j < adjlist[i].size(); j++) {
				if (len(pos(i), pos(adjlist[i][j])) == 0) {
					// if the weight has not been calculated 
					len(pos(i), pos(adjlist[i][j])) = sqrt(
						(V(i, 0) - V(adjlist[i][j], 0))*(V(i, 0) - V(adjlist[i][j], 0)) +
						(V(i, 1) - V(adjlist[i][j], 1))*(V(i, 1) - V(adjlist[i][j], 1)) +
						(V(i, 2) - V(adjlist[i][j], 2))*(V(i, 2) - V(adjlist[i][j], 2))
					);
					len(pos(adjlist[i][j]), pos(i)) = len(pos(i), pos(adjlist[i][j]));
				}
				sum_wi(pos(i)) += len(pos(i), pos(adjlist[i][j]));
			}
		}
	}

	// construct A 
	for (int i = 0; i < A.rows(); i++) {
		for (int j = 0; j < A.cols(); j++) {
			// construct A[i][j]
			if (j == i) {
				// A[i][j] = 1 + ..
				A(i, j) = 1;
				for (int k = 0; k < len.cols(); k++) {
					if (sum_wi(k) != 0) {
						// due to construction of sum_wi, sum_wi(onering) = 0				
						A(i, j) += (1 / sum_wi(pos(num_new + i)))*((len(pos(num_new + i), k)*len(pos(num_new + i), k)) / sum_wi(k));
					}
				}
			}
			else {
				// A[i][j] = -2 *.. + ..
				A(i, j) = -2 * len(pos(num_new + i), pos(num_new + j)) / sum_wi(pos(num_new + j));
				for (int k = 0; k < len.cols(); k++) {
					if (sum_wi(k) != 0) {
						// due to construction of sum_wi, sum_wi(onering) = 0
						A(i, j) += (1 / sum_wi(pos(num_new + i)))*((len(pos(num_new + i), k)*len(k, pos(num_new + j))) / sum_wi(k));
					}
				}
			}
		}
	}

	//double tmo;
	// construct B
	for (int i = 0; i < B.rows(); i++) {
		// initialize B(x,y,z)
		B(i, 0) = 0;
		B(i, 1) = 0;
		B(i, 2) = 0;
		for (int j = 0; j < len.cols(); j++) {
			if (posmark(pos_back(j)) == 1 || posmark(pos_back(j)) == 2) {
				// only calc with boudnary and onering vertex
				for (int k = 0; k < len.cols(); k++) {
					if (sum_wi(k) != 0) {
						B(i, 0) += -(1 / sum_wi(pos(num_new + i)))*(len(pos(num_new + i), k)*len(k, j)*V(pos_back(j), 0)) / sum_wi(k);
						B(i, 1) += -(1 / sum_wi(pos(num_new + i)))*(len(pos(num_new + i), k)*len(k, j)*V(pos_back(j), 1)) / sum_wi(k);
						B(i, 2) += -(1 / sum_wi(pos(num_new + i)))*(len(pos(num_new + i), k)*len(k, j)*V(pos_back(j), 2)) / sum_wi(k);	
						//tmo = -(1 / sum_wi(pos(num_new + i)))*(len(pos(num_new + i), k)*len(k, j)*V(pos_back(j), 0)) / sum_wi(k);
						//if (i == B.rows() - 1 && tmo != 0) {
						//	std::cout << "#j = " << pos_back(j);
						//	std::cout << "#k = " << pos_back(k);
						//	std::cout << "-----1 / " << sum_wi(pos(num_new + i)) << " * " << len(pos(num_new + i), k) << " * " <<len(k,j)<<" * "<< V(pos_back(j), 0) << std::endl;
						//}
					}
				}
				B(i, 0) += (2 / sum_wi(pos(num_new + i)))*len(pos(num_new + i), j)*V(pos_back(j), 0);
				B(i, 1) += (2 / sum_wi(pos(num_new + i)))*len(pos(num_new + i), j)*V(pos_back(j), 1);
				B(i, 2) += (2 / sum_wi(pos(num_new + i)))*len(pos(num_new + i), j)*V(pos_back(j), 2);

			}
		}
	}
	X = A.colPivHouseholderQr().solve(B);
	int kk = 0;
	for (int i = num_new; i < V.rows(); i++) {
		std::cout << "v xyz " << V.row(i) << std::endl;
		V(i, 0) += X(kk, 0);

		V(i, 1) += X(kk, 1);
		V(i, 2) += X(kk, 2);
		std::cout << "v xyz " << V.row(i) << std::endl;

		kk++;
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
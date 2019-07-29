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

	// initialize len and sum_wi  
	len.resize(num_total, num_total); // len(i, j) = len||vi - vj||
	len.setZero();
	sum_wi.resize(num_total);
	sum_wi.setZero();
	pos_back.resize(num_total);

	// initialize AX=B
	A.resize(count_new, count_new);
	A.setIdentity();
	X.resize(count_new, 3);
	B.resize(count_new, 3);
	B.setZero();


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
					// [NOTE] the scale-dependent laplacian operator should be 1/||vi-vj|| which is wrong in Nan Qiu's paper(in this paper is ||vi - vj||)
					// in Liepa's paper is 1/||vi - vj||
					len(pos(i), pos(adjlist[i][j])) = 1/(sqrt(
						(V(i, 0) - V(adjlist[i][j], 0))*(V(i, 0) - V(adjlist[i][j], 0)) +
						(V(i, 1) - V(adjlist[i][j], 1))*(V(i, 1) - V(adjlist[i][j], 1)) +
						(V(i, 2) - V(adjlist[i][j], 2))*(V(i, 2) - V(adjlist[i][j], 2))
					));
					len(pos(adjlist[i][j]), pos(i)) = len(pos(i), pos(adjlist[i][j]));
				}
				sum_wi(pos(i)) += len(pos(i), pos(adjlist[i][j]));
			}
		}
	}
	//std::ofstream ff1,ff2;
	//ff1.open("log.txt");
	//ff2.open("Aconstruct.txt");
	// construct A 
	for (int i = 0; i < A.rows(); i++) {
		//ff1 << "construct line ====== " << i << " coresponding to original vertex " << num_new + i << std::endl;
		for (int j = 0; j < A.cols(); j++) {
			// construct A[i][j]
			if (j == i) {
				// A[i][j] = 1 + ..
				for (int k = 0; k < sum_wi.cols(); k++) {
					if ((sum_wi(k) != 0) && (len(pos(num_new + i), k))) {
						// due to construction of sum_wi, sum_wi(onering) = 0				
						A(i, j) += (1 / sum_wi(pos(num_new + i)))*((len(pos(num_new + i), k)*len(pos(num_new + i), k)) / sum_wi(k));
						//ff1 << " k = " << k << " coresponding to original vertex " << pos_back(k) << std::endl;
						//ff1 << "calc A( " << i << " , " << j << " ) = " << (1 / sum_wi(pos(num_new + i))) << " * " << len(pos(num_new + i), k)
						//	<< " * " << len(pos(num_new + i), k) << " / " << sum_wi(k) << std::endl;
					}
				}
			}
			else {
				//ff1 << "[initialization] i, j corresponding to " << num_new + i << " , " << num_new + j << std::endl;
				//ff1 << std::endl;

				if (len(pos(num_new + i), pos(num_new + j)) != 0) {
					A(i, j) = -2 * len(pos(num_new + i), pos(num_new + j)) / sum_wi(pos(num_new + i));
					//ff1 << " calc A( " << i << " , " << j << " )= " << " -2 * " << len(pos(num_new + i), pos(num_new + j)) << " / " << sum_wi(pos(num_new + i)) << std::endl;
				}
				// A[i][j] = -2 *.. + ..

				for (int k = 0; k < len.cols(); k++) { 
					if ((sum_wi(k) != 0) && (len(pos(num_new + i), k) != 0) && (len(k, pos(num_new + j)) != 0)) {

						// due to construction of sum_wi, sum_wi(onering) = 0
						A(i, j) += (1 / sum_wi(pos(num_new + i)))*((len(pos(num_new + i), k)*len(k, pos(num_new + j))) / sum_wi(k));
						//ff1 << " k = " << k << " coresponding to original vertex " << pos_back(k) << std::endl;
						//ff1 << "calc A( " << i << " , " << j << " ) += " << (1 / sum_wi(pos(num_new + i))) << " * " << len(pos(num_new + i), k)
						//	<< " * " << len(k, pos(num_new + j)) << " / " << sum_wi(k) << std::endl;
					}
				}
			}
		}
	}

	//ff1.close();
	//ff2.close();
	//double tmo;

	//std::ofstream ff3, ff4;
	////ff3.open("log.txt");
	//ff4.open("Bconstruct.txt");
	// construct B
	for (int i = 0; i < B.rows(); i++) {
		//ff3 << "construct line ====== " << i << " coresponding to original vertex " << num_new + i << std::endl;
		//ff3 << std::endl;
		for (int j = 0; j < num_total-count_new; j++) {

			// only calc with boudnary and onering vertex
			// B left = 2 * 1/sum(line i)*w(vi, vj)*V(j)
			if (len(pos(num_new + i), j) != 0) {
				//ff3 << "calc [LEFT] vertex j =  " << j << " coresponding to original vertex " << pos_back(j) << std::endl;

				B(i, 0) += 2 * len(pos(num_new + i), j)*V(pos_back(j), 0) / sum_wi(pos(num_new + i));
				B(i, 1) += 2 * len(pos(num_new + i), j)*V(pos_back(j), 1) / sum_wi(pos(num_new + i));
				B(i, 2) += 2 * len(pos(num_new + i), j)*V(pos_back(j), 2) / sum_wi(pos(num_new + i));
				//ff3 << "calc [LEFT] x<===>B(i,0) += " << " 2 * " << len(pos(num_new + i), j) << " * " << V(pos_back(j), 0) << " / " << sum_wi(pos(num_new + i)) << std::endl;
				//ff3 << "calc [LEFT] y<===>B(i,1) += " << " 2 * " << len(pos(num_new + i), j) << " * " << V(pos_back(j), 1) << " / " << sum_wi(pos(num_new + i)) << std::endl;
				//ff3 << "calc [LEFT] z<===>B(i,2) += " << " 2 * " << len(pos(num_new + i), j) << " * " << V(pos_back(j), 2) << " / " << sum_wi(pos(num_new + i)) << std::endl;
				//ff3 << std::endl;
			}
			// B right = -1/sum(line i) * ..
			for (int k = 0; k < sum_wi.cols(); k++) {
				if ((sum_wi(k) != 0) && (len(pos(num_new + i), k) != 0) && (len(k, j) != 0)) {
					//ff3 << "calc [RIGHT] vertex j =  " << j << " coresponding to original vertex " << pos_back(j) << std::endl;
					//ff3 << "calc [RIGHT] k = " << k << " corresponding to original vertex " << pos_back(k) << std::endl;
					//ff3 << " calc w( " << num_new + i << " , " << pos_back(k) << " ) * w( " << pos_back(k) << " , " << pos_back(j) << std::endl;
					B(i, 0) += -(1 / sum_wi(pos(num_new + i)))*(len(pos(num_new + i), k)*len(k, j)*V(pos_back(j), 0)) / sum_wi(k);
					B(i, 1) += -(1 / sum_wi(pos(num_new + i)))*(len(pos(num_new + i), k)*len(k, j)*V(pos_back(j), 1)) / sum_wi(k);
					B(i, 2) += -(1 / sum_wi(pos(num_new + i)))*(len(pos(num_new + i), k)*len(k, j)*V(pos_back(j), 2)) / sum_wi(k);
					//ff3 << "calc x<===>B(i,0) += " << " ( -1 / " << sum_wi(pos(num_new + i)) << " ) * " <<
					//	len(pos(num_new + i), k) << " * " << len(k, j) << " * " << V(pos_back(j), 0) << " / " << sum_wi(k) << std::endl;
					//ff3 << "calc y<===>B(i,1) += " << " ( -1 / " << sum_wi(pos(num_new + i)) << " ) * " <<
					//	len(pos(num_new + i), k) << " * " << len(k, j) << " * " << V(pos_back(j), 1) << " / " << sum_wi(k) << std::endl;
					//ff3 << "calc z<===>B(i,2) += " << " ( -1 / " << sum_wi(pos(num_new + i)) << " ) * " <<
					//	len(pos(num_new + i), k) << " * " << len(k, j) << " * " << V(pos_back(j), 2) << " / " << sum_wi(k) << std::endl;
				}
			}
		}
		//ff3 << std::endl;
	}

	//ff3.close();
	//ff4 << B << std::endl;
	//ff4.close();

	X = A.colPivHouseholderQr().solve(B);
	//std::ofstream ff5;
	//ff5.open("X.txt");
	//ff5 << X << std::endl;
	//ff5.close();
	int kk = 0;
	for (int i = num_new; i < V.rows(); i++) {
		V.row(i) = X.row(kk);
		std::cout << X.row(kk) << std::endl;
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
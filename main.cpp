/*
	Created by Yuhan Ping from HKU at 2019-06-04.
	The project is about tooth completion.
	Under development.
	Lateset updated on 2019-06-04.
*/


#include <iostream>
#include <igl/readOFF.h>
#include "Header/io.h"
#include <igl/boundary_loop.h>
#include <igl/fit_plane.h>

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/ml.hpp>

using namespace Eigen;
using namespace std;
using namespace cv;
using namespace cv::ml;

// ---------------pre-define parameters------------------
MatrixXd V_plane, V_total, V1_boundary_L, V1_boundary_R, V2_boundary_L, V2_boundary_R;

RowVectorXd L1, L2;
RowVector3d N, C; 


string first_mesh = "F:/StudyMaterials/HKU/RA/libigl-example-project/32.off";
string second_mesh = "F:/StudyMaterials/HKU/RA/libigl-example-project/33.off";
/*
	N: 1 * 3 vector, the normal of fitted plane
	C: 1 * 3 vector, a point lies in the plane
*/
//--------------------------------------------------------

int main(int argc, char *argv[])
{
	// load the mesh, currently load one mesh for test, load two mesh in the future for dealing with neighbor intersection
	load_mesh(first_mesh,V1,F1);

	// identify the boundary
	//igl::boundary_loop(F1, L1);
	//cout << "boundary vertices idx " << L1 << endl;
	
	// create interactive straight line

	// TODO: triangulate the hole 

	// TODO: refinement

	// TODO: fairing with the surface conditions

	// Plot the mesh
	viewer_display(V1, F1);
	
	getchar();
	return 0;
	// test opencv svm in 2D 
	// setup training data
	//int labels[4] = { 1, -1, -1, 1 };
	//float trainingData[4][2] = { {501,10}, {255,10},{501,255},{10,501} };
	//Mat trainingDataMat(4, 2, CV_32F, trainingData);
	//Mat labelsMat(4, 1, CV_32SC1, labels);

	//// train the svm
	//Ptr<SVM> svm = SVM::create();
	//svm->setType(SVM::C_SVC);
	//svm->setKernel(SVM::LINEAR);
	//svm->setTermCriteria(TermCriteria(TermCriteria::MAX_ITER, 100, 1e-6));
	//svm->train(trainingDataMat, ROW_SAMPLE, labelsMat);

	//Mat sampleMat = (Mat_<float>(1, 2) << 1, 2);
	//// data for visual representation
	//int width = 512;
	//int height = 512;
	//Mat image = Mat::zeros(height, width, CV_8UC3);

	//// show the decision regions given by the svm
	//Vec3b green(0, 255, 0);
	//Vec3b blue(0, 0, 255);
	//for (int i = 0; i < image.rows; i++) {
	//	for (int j = 0; i < image.cols; j++) {
	//		Mat sampleMat = (Mat_<float>(1, 2) << j, i);
	//		float response = svm->predict(sampleMat);
	//	}
	//}
	


	//// load the mesh
	//igl::readOFF("F:/StudyMaterials/HKU/RA/libigl-example-project/32.off", V1, F1);
	//igl::readOFF("F:/StudyMaterials/HKU/RA/libigl-example-project/33.off", V2, F2);

	////cout << "v1 rows v2 rows v1 cols v2 cols " << V1.rows() << " " << V2.rows() << " " << V1.cols() << " " << V2.cols() << endl;
	//
	//// test cube.off: 
	////igl::readOFF("F:/StudyMaterials/HKU/RA/libigl-example-project/cube.off", V1, F1);
	//// get boundary vertices and faces
	//igl::boundary_loop(F1, L1);
	//igl::boundary_loop(F2, L2);
	//cout << "L1 L2 " << L1 <<endl<<L2<< endl; // L stores the index of boundary vertex with order
	////cout << "L1 rows cols" << L1.rows() <<" "<< L1.cols() << endl;

	//double xmid = 0.5 * (V1.col(0).maxCoeff() + V1.col(0).minCoeff());
	//cout << "xmax = " << V1.col(0).maxCoeff() << endl;
	//cout << "xmin = " << V1.col(0).minCoeff() << endl;
	//cout << "xmid = " << xmid << endl;

	////for (int i = 0; i < L1.rows(); i++) {
	////	if(V1(L1(0,i))>0.5*)
	////}

	//// TO BE DISCUSSED: smooth the boundary and remove unnecessary faces?

	//// TODO: get a plane between to teeth boundarys 
	//// use 
	// igl::fit_plane(V_total, N, C);
	// cout << "normal " << N << endl;
	// cout << "one point " << C << endl;

	// V_plane.resize(4, 3);
	// cout << "v_plane size " << V_plane.rows() << " " << V_plane.cols() << endl;

	// // x = -1/A * (B(y -y0)+C(z-z0)) + x0
	// double x0 = C(0, 0);
	// cout << "x0 " << x0 << endl;
	// double ymin = V1.col(1).minCoeff() - 1;
	// double ymax = V1.col(1).maxCoeff() + 1;
	// double zmin = V1.col(2).minCoeff() - 1;
	// double zmax = V1.col(2).maxCoeff() + 1;
	// cout << "ymin ymax zmin zmax " << ymin << " " << ymax << " " << zmin << " " << zmax << endl;

	// double x_p0 = x0 - (N(0, 1)*(ymin - C(0, 1)) + N(0, 2)*(zmin - C(0, 2))) / N(0, 0);
	// double x_p1 = x0 - (N(0, 1)*(ymin - C(0, 1)) + N(0, 2)*(zmax - C(0, 2))) / N(0, 0);
	// double x_p2 = x0 - (N(0, 1)*(ymax - C(0, 1)) + N(0, 2)*(zmax - C(0, 2))) / N(0, 0);
	// double x_p3 = x0 - (N(0, 1)*(ymax - C(0, 1)) + N(0, 2)*(zmin - C(0, 2))) / N(0, 0);



	// V_plane <<
	//	 x_p0, ymin, zmin,
	//	 x_p1, ymin, zmax,
	//	 x_p2, ymax, zmax,
	//	 x_p3, ymax, zmin;

	// cout << "v plane " << V_plane << endl;

	// //V_plane(0, 0) = C.col(0) + ((N.col(1) * (ymin - C.col(1)) + N.col(2) *(zmin - C.col(2))) / N.col(0)) * -1; //p1
	// //V_plane(1, 0) = C.col(0) + ((N.col(1) * (ymax - C.col(1)) + N.col(2) *(zmin - C.col(2))) / N.col(0)) * -1; //p2
	// //V_plane(2, 0) = C.col(0) + ((N.col(1) * (ymin - C.col(1)) + N.col(2) *(zmax - C.col(2))) / N.col(0)) * -1; //p3
	// //V_plane(3, 0) = C.col(0) + ((N.col(1) * (ymax - C.col(1)) + N.col(2) *(zmax - C.col(2))) / N.col(0)) * -1; //p4
	// 
	// //V_plane.resize(4, 3);
	// //double 

	//// TODO: triangulate the hole 

	//// TODO: refinement

	//// TODO: fairing with the surface conditions

	//// Plot the mesh

	//igl::opengl::glfw::Viewer viewer;
	//igl::opengl::glfw::imgui::ImGuiMenu menu;
	//viewer.plugins.push_back(&menu);
	//viewer.data().set_mesh(V1, F1);
	//viewer.data().add_points(C, RowVector3d(0.0, 255.0, 0.0));
	//viewer.data().add_points(V_plane, RowVector3d(255.0, 0.0, 0.0));
	//viewer.data().add_edges(V_plane.row(0), V_plane.row(1), RowVector3d(255.0, 0.0, 0.0));
	//viewer.data().add_edges(V_plane.row(1), V_plane.row(2), RowVector3d(255.0, 0.0, 0.0));
	//viewer.data().add_edges(V_plane.row(2), V_plane.row(3), RowVector3d(255.0, 0.0, 0.0));
	//viewer.data().add_edges(V_plane.row(3), V_plane.row(0), RowVector3d(255.0, 0.0, 0.0));
	//viewer.append_mesh();
	//viewer.data().set_mesh(V2, F2);

 // viewer.launch();
}

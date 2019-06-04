#include "Header/io.h"

Eigen::Vector3f select_v1, select_v2;
Eigen::MatrixXd V1;
Eigen::MatrixXi F1;

bool mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
	//, Eigen::Vector3f &v1, Eigen::Vector3f &v2
	// cast a ray in the view direction starting from the mouse position to get the real 3D position on the mesh
	int fid;
	Eigen::Vector3f vid;
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;

	std::cout << "current mouse location-----[" << x << ", " << y << "]" << std::endl;

	if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view, viewer.core.proj, viewer.core.viewport, V1, F1, fid, vid)) {
		std::cout << "select face " << fid << std::endl;
		std::cout << "select vertex " << vid << std::endl;
	}


	return true;
}


bool viewer_display(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);
	viewer.callback_mouse_down = &mouse_down;
	viewer.data().set_mesh(V, F);
	viewer.launch();

	return true;
}; 

bool load_mesh(std::string &mesh_dir, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	igl::readOFF(mesh_dir, V, F);
	std::cout << "Load mesh successfully----[" << mesh_dir << "]" << std::endl;
	return true;
}
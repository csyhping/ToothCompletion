#include "Header/io.h"
#include "Header/boundary_construction.h"

Eigen::RowVector3d select_v1, select_v2;
Eigen::MatrixXd V1;
Eigen::MatrixXi F1;
int select_count = 0; 




bool mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier) {
	
	// cast a ray in the view direction starting from the mouse position to get the real 3D position on the mesh
	int fid; // selected face id 
	Eigen::Vector3f vid; // barycentric coordinate of selected face, vid = [a1, a2, a3] which represents the point p = a1*v1 + a2*v2 + a3*v3
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;

	std::cout << "current mouse location-----[" << x << ", " << y << "]" << std::endl; 

	if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view, viewer.core.proj, viewer.core.viewport, V1, F1, fid, vid)) {
		long c; // which vertex has the max value, a1/a2/a3? which means c = 0 or 1 or 2 of the select face
		vid.maxCoeff(&c);
		switch (select_count)
		{
		case 0:
			// locate the nearest vertex as selected v1
			select_v1 = V1.row(F1(fid, c));
			std::cout << "the maxcoeff of vid " << c << std::endl;
			std::cout << "[INFO]----First vertex selected: [" << select_v1 << "]" << std::endl;
			
			// visualize v1
			viewer.data().add_points(select_v1,Eigen::RowVector3d(255,0,0));
			select_count++;
			break;
		case 1:
			// locate the nearest vertex as selected v2
			select_v2 = V1.row(F1(fid, c));
			std::cout << "the maxcoeff of vid " << c << std::endl;
			std::cout << "[INFO]----Second vertex selected: [" << select_v2 << "]" << std::endl;

			// visualize v2
			viewer.data().add_points(select_v2, Eigen::RowVector3d(0, 255, 0));
			select_count++;
			break;
		default:
			std::cout << "[INFO]----Already selected two vertex" << std::endl;
			break;

		}
		// only if the mouse is on the mesh, the mouse is able to select vertex
		return true;
	}
	
	// if the mouse is not on mesh, the mouse is able to rotate the mesh as normal
	return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
	if (key == ' ') {
		// presee [space] to go further
		// create line
		viewer.data().add_edges(select_v1, select_v2, Eigen::RowVector3d(0, 0, 255));
		return true;
	}
	return false;
}

bool viewer_display(Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);
	viewer.callback_mouse_down = &mouse_down;
	viewer.callback_key_down = &key_down;
	viewer.data().set_mesh(V, F);
	viewer.data().set_colors(Color_per_vertex);
	viewer.launch();

	return true;
}; 

bool load_mesh(std::string &mesh_dir, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	igl::readOFF(mesh_dir, V, F);
	std::cout << "Load mesh successfully----[" << mesh_dir << "]" << std::endl;
	return true;
}
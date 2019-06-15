#include "Header/io.h"
#include "Header/boundary_construction.h"
#include "Header/plane.h"

Eigen::RowVector3d select_v1, select_v2, select_v3, select_v4;
Eigen::MatrixXd V1, New_vertex_on_line_R, New_vertex_on_line_L;
Eigen::MatrixXd Projected_vertex_R, Projected_vertex_L; // projected vertex on 2D plane, #row = hole_vertex, #col = 2 (x,y)
Eigen::MatrixXd Vertex_on_xy_R, Vertex_on_xy_L; // rotated vertex on the xy plane
Eigen::MatrixXd Hole_vertex_R, Hole_vertex_L;// hole_boundary vertices, including select_v1 + select_v2 + new_v_on_line + orginal_boundary_v_above_v1&v2
Eigen::MatrixXi F1;
Eigen::RowVector3d NR, NL, CR, CL; // NR & NL: the normal of the right & left plane, CR & CL: one point on the right & left plane
int select_count_L = 0; 
int select_count_R = 0;
int count_L = 0;
int count_R = 0;
int idx_v1, idx_v2, idx_v3, idx_v4;




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
		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			// left click to select two vertex

			switch (select_count_R)
			{
			case 0:
				// locate the nearest vertex as selected v1
				select_v1 = V1.row(F1(fid, c));
				idx_v1 = F1(fid, c);
				std::cout << "the maxcoeff of vid " << c << std::endl;
				std::cout << "[INFO]----First vertex selected: [" << select_v1 << "]" << " ------V_" << idx_v1 << std::endl;


				// visualize v1
				viewer.data().add_points(select_v1, Eigen::RowVector3d(255, 0, 0));
				select_count_R++;
				break;
			case 1:
				// locate the nearest vertex as selected v2
				select_v2 = V1.row(F1(fid, c));
				idx_v2 = F1(fid, c);
				std::cout << "the maxcoeff of vid " << c << std::endl;
				std::cout << "[INFO]----Second vertex selected: [" << select_v2 << "]" << " ------V_" << idx_v2 << std::endl;

				// visualize v2
				viewer.data().add_points(select_v2, Eigen::RowVector3d(0, 255, 0));
				select_count_R++;
				break;
			default:
				std::cout << "[INFO]----Already selected two vertex" << std::endl;
				break;

			}
			// only if the mouse is on the mesh, the mouse is able to select vertex
			return true;
		}
		else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
			// right click to select two vertex

			switch (select_count_L)
			{
			case 0:
				// locate the nearest vertex as selected v1
				select_v3 = V1.row(F1(fid, c));
				idx_v3 = F1(fid, c);
				std::cout << "the maxcoeff of vid " << c << std::endl;
				std::cout << "[INFO]----First vertex selected: [" << select_v3 << "]" << " ------V_" << idx_v3 << std::endl;


				// visualize v1
				viewer.data().add_points(select_v3, Eigen::RowVector3d(255, 0, 0));
				select_count_L++;
				break;
			case 1:
				// locate the nearest vertex as selected v2
				select_v4 = V1.row(F1(fid, c));
				idx_v4 = F1(fid, c);
				std::cout << "the maxcoeff of vid " << c << std::endl;
				std::cout << "[INFO]----Second vertex selected: [" << select_v4 << "]" << " ------V_" << idx_v4 << std::endl;

				// visualize v2
				viewer.data().add_points(select_v4, Eigen::RowVector3d(0, 255, 0));
				select_count_L++;
				break;
			default:
				std::cout << "[INFO]----Already selected two vertex" << std::endl;
				break;

			}
			// only if the mouse is on the mesh, the mouse is able to select vertex
			return true;
		}

	}
	
	// if the mouse is not on mesh, the mouse is able to rotate the mesh as normal
	return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
	// [NOTE]: can add more functions here to go on process
	if (key == ' ') {
		// presee [space] to go further
		// create line
		viewer.data().add_edges(select_v1, select_v2, Eigen::RowVector3d(0, 0, 255)); // between right two selected_v
		viewer.data().add_edges(select_v3, select_v4, Eigen::RowVector3d(0, 0, 255)); // between left two selected_v

		
		// get new vertex coordinates
		create_vertex_on_line(select_v1, select_v2, New_vertex_on_line_R, count_R); // for right edge
		create_vertex_on_line(select_v3, select_v4, New_vertex_on_line_L, count_L); // for left edge

		// visualize new vertex on line
		viewer.data().add_points(New_vertex_on_line_R, Eigen::RowVector3d(0, 255, 255)); // for right edge
		viewer.data().add_points(New_vertex_on_line_L, Eigen::RowVector3d(0, 255, 255)); // for left edge

		// get new hole boundary
		get_hole_boundary(V1, F1, select_v1, select_v2, New_vertex_on_line_R, idx_v1, idx_v2, count_R,Hole_vertex_R); // for right hole
		get_hole_boundary(V1, F1, select_v3, select_v4, New_vertex_on_line_L, idx_v3, idx_v4, count_L,Hole_vertex_L); // for right hole

		// fit plane
		get_plane(Hole_vertex_R, NR, CR); // for right side
		get_plane(Hole_vertex_L, NL, CL); // for left side

		// project hole vertex to plane
		project_hole_vertex_to_plane(Hole_vertex_R, NR, CR, Projected_vertex_R);
		project_hole_vertex_to_plane(Hole_vertex_L, NL, CL, Projected_vertex_L);

		// rotate the plane to a xy plane
		rotate_to_xy_plane(NR, Projected_vertex_R, Vertex_on_xy_R);

		// test the projected vertex
		viewer.data().add_points(Projected_vertex_R, Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_points(Projected_vertex_L, Eigen::RowVector3d(255, 255, 0));

		// test the plane 
		viewer.data().add_points(CR, Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_points(CL, Eigen::RowVector3d(255, 255, 0));


		Eigen::MatrixXd test_R(4, 3);
		Eigen::MatrixXd test_L(4, 3);
		double ymaxR = CR(0, 1) + 2;
		double yminR = CR(0, 1) - 2;
		double zminR = CR(0, 2) - 2;
		double zmaxR = CR(0, 2) + 2;
		double ymaxL = CL(0, 1) + 2;
		double yminL = CL(0, 1) - 2;
		double zminL = CL(0, 2) - 2;
		double zmaxL = CL(0, 2) + 2;

		test_R.row(0) <<
			CR(0, 0) - (NR(0, 1)*(yminR - CR(0, 1)) + NR(0, 2)*(zminR - CR(0, 2))) / NR(0, 0), yminR, zminR;
		test_R.row(1) <<
			CR(0, 0) - (NR(0, 1)*(ymaxR - CR(0, 1)) + NR(0, 2)*(zminR - CR(0, 2))) / NR(0, 0), ymaxR, zminR;
		test_R.row(2) <<
			CR(0, 0) - (NR(0, 1)*(yminR - CR(0, 1)) + NR(0, 2)*(zmaxR - CR(0, 2))) / NR(0, 0), yminR, zmaxR;
		test_R.row(3) <<
			CR(0, 0) - (NR(0, 1)*(ymaxR - CR(0, 1)) + NR(0, 2)*(zmaxR - CR(0, 2))) / NR(0, 0), ymaxR, zmaxR;
		test_L.row(0) <<
			CL(0, 0) - (NL(0, 1)*(yminL - CL(0, 1)) + NL(0, 2)*(zminL - CL(0, 2))) / NL(0, 0), yminL, zminL;
		test_L.row(1) <<
			CL(0, 0) - (NL(0, 1)*(ymaxL - CL(0, 1)) + NL(0, 2)*(zminL - CL(0, 2))) / NL(0, 0), ymaxL, zminL;
		test_L.row(2) <<
			CL(0, 0) - (NL(0, 1)*(yminL - CL(0, 1)) + NL(0, 2)*(zmaxL - CL(0, 2))) / NL(0, 0), yminL, zmaxL;
		test_L.row(3) <<
			CL(0, 0) - (NL(0, 1)*(ymaxL - CL(0, 1)) + NL(0, 2)*(zmaxL - CL(0, 2))) / NL(0, 0), ymaxL, zmaxL;

		viewer.data().add_points(test_R, Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CR, test_R.row(0), Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CR, test_R.row(1), Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CR, test_R.row(2), Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CR, test_R.row(3), Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_points(test_L, Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CL, test_L.row(0), Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CL, test_L.row(1), Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CL, test_L.row(2), Eigen::RowVector3d(255, 255, 0));
		viewer.data().add_edges(CL, test_L.row(3), Eigen::RowVector3d(255, 255, 0));




		viewer.data().set_colors(Color_per_vertex);

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
	// visualization test, color matrix initialization
	Color_per_vertex = Eigen::MatrixXd::Constant(V.rows(), 3, 1);
	viewer.data().set_colors(Color_per_vertex);
	viewer.launch();

	return true;
}; 

bool load_mesh(std::string &mesh_dir, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	igl::readOFF(mesh_dir, V, F);
	std::cout << "Load mesh successfully----[" << mesh_dir << "]" << std::endl;
	get_pos_boundary(F);
	return true;
}
#include "half_edge_simp/half_edge_mesh.h"

#include <iostream>
#include <gflags/gflags.h>
#include <igl/readPLY.h>

DEFINE_string(ply_fn, "000202.ply", "ply file large faces");
DEFINE_int32(target_count, 30000, "ply file large faces");

#pragma comment(lib, "igl.lib")

int main(int argc, char** argv)
{
	std::string dat_dir = argv[1];
	std::string ply_fn = dat_dir + "/" + FLAGS_ply_fn;
	std::string result_fn = dat_dir + "/result.obj";
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::readPLY(ply_fn.c_str(), V, F);

	int target_count = FLAGS_target_count;
	half_edge::Mesh mesh;
	mesh.simplifyMesh(V, F, target_count, result_fn);

	V = mesh.getSimplifiedV();
	F = mesh.getSimplifiedF();

	printf("done.\n");
	return 0;
}
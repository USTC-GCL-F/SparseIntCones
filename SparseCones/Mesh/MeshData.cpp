#include "MeshData.h"
#include <fstream>

MeshData::MeshData(Mesh& mesh) :mesh(mesh)
{

}

void MeshData::VertArea()
{
	A_vec.setConstant(mesh.n_vertices(), 0);
	double fA;
	for (const FH& f_h : mesh.faces())
	{
		fA = mesh.calc_face_area(f_h) / 3;
		for (const VH& fv_h : mesh.fv_range(f_h))
		{
			A_vec[fv_h.idx()] += fA;
		}
	}

	A_vec = A_vec / A_vec.sum();

	/*
	std::ofstream _out("A_arm1.txt");
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		_out << A_vec[i] << std::endl;
	}
	_out.close();*/

}

void MeshData::FaceAngle()
{
	f_angle.setConstant(mesh.n_faces() * 3, 0);
	for (const FH& f_h : mesh.faces())
	{
		std::vector<OpenMesh::Vec3d> he_v(3);
		int cout = 0;
		for (const HEH& he_h : mesh.fh_range(f_h))
		{
			he_v[cout] = mesh.calc_edge_vector(he_h).normalized();
			cout++;
		}

		double cos_0 = -dot(he_v[1], he_v[2]);
		double cos_1 = -dot(he_v[2], he_v[0]);
		double cos_2 = -dot(he_v[0], he_v[1]);

		f_angle[3 * f_h.idx() + 0] = acos(cos_0);
		f_angle[3 * f_h.idx() + 1] = acos(cos_1);
		f_angle[3 * f_h.idx() + 2] = acos(cos_2);
	}
}

void MeshData::VertGauss()
{
	K_vec.setConstant(mesh.n_vertices(), 2 * M_PI);

	for (const VH& v_h : mesh.vertices())
	{
		if (mesh.is_boundary(v_h)) {
			K_vec[v_h.idx()] = M_PI;
		}
	}

	for (const FH& f_h : mesh.faces())
	{
		std::vector<int> v_ids(3);
		int cout = 0;
		for (const HEH& he_h : mesh.fh_range(f_h))
		{
			v_ids[cout] = mesh.to_vertex_handle(mesh.next_halfedge_handle(he_h)).idx();
			cout++;
		}

		K_vec[v_ids[0]] -= f_angle[3 * f_h.idx() + 0];
		K_vec[v_ids[1]] -= f_angle[3 * f_h.idx() + 1];
		K_vec[v_ids[2]] -= f_angle[3 * f_h.idx() + 2];
	}

	/*
	std::ofstream _out("K_ori_arm1.txt");
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		_out << K_vec[i] << std::endl;
	}
	_out.close();*/

}


void MeshData::MatrixL()
{
	double cot_weight;
	VH to_v, from_v;
	std::vector<T> trip;
	trip.reserve(12 * mesh.n_faces());
	for (const FH& f_h : mesh.faces())
	{
		int cout = 0;
		for (const HEH& fh_h : mesh.fh_range(f_h))
		{
			to_v = mesh.to_vertex_handle(fh_h);
			from_v = mesh.from_vertex_handle(fh_h);

			cot_weight = 0.5 / tan(f_angle[3 * f_h.idx() + cout]);
			cout++;

			trip.emplace_back(to_v.idx(), from_v.idx(), -cot_weight);
			trip.emplace_back(from_v.idx(), to_v.idx(), -cot_weight);
			trip.emplace_back(from_v.idx(), from_v.idx(), cot_weight);
			trip.emplace_back(to_v.idx(), to_v.idx(), cot_weight);
		}
	}

	L_matrix.resize(mesh.n_vertices(), mesh.n_vertices());
	L_matrix.setFromTriplets(trip.begin(), trip.end());
	/*
	std::ofstream _out("L_arm1.txt");
	for (int k = 0; k < L_matrix.outerSize(); ++k)
	{
		for (SparseMatrix<double>::InnerIterator it(L_matrix, k); it; ++it)
		{
			_out << it.row() + 1 << " " << it.col() + 1 << " " << it.value() << std::endl;

		}
	}*/

}

void MeshData::MatrixP()
{
	P.resize(mesh.n_vertices(), mesh.n_vertices() - 1);
	std::vector<T> element;
	for (int i = 0; i < mesh.n_vertices() - 1; i++)
	{
		element.emplace_back(i, i, 1);
	}

	P.setFromTriplets(element.begin(), element.end());
	PT = P.transpose();
}

void MeshData::Init()
{
	VertArea();
	FaceAngle();
	VertGauss();
	MatrixL();
	MatrixP();
}
#include "QEMSimplification.h"
#include <iostream>

QEMSimplification::QEMSimplification(const Mesh & m)
	:mesh(m)
{
	mesh.add_property(vindex);
	for (const auto & vh : mesh.vertices())
	{
		mesh.property(vindex, vh) = vh.idx();
	}
	mesh.add_property(QuaMat_f);
	mesh.add_property(QuaMat_v);
	mesh.add_property(QuaMat_e);
	mesh.add_property(area_f);
	mesh.add_property(new_point);
	mesh.add_property(is_cal_e);
	mesh.add_property(is_cal_he);
	for (const auto & eh : mesh.edges())
	{
		mesh.property(is_cal_e, eh) = false;
		mesh.property(is_cal_he, mesh.halfedge_handle(eh, 0)) = false;
		mesh.property(is_cal_he, mesh.halfedge_handle(eh, 1)) = false;
	}
}

QEMSimplification::~QEMSimplification(void)
{
}

void QEMSimplification::Simplify(int tar_num_v, double threshold)
{
	InitialEdgeCost();
	index_collapse_v_to_v.clear();
	index_collapse_v_to_v.reserve(mesh.n_vertices());
	auto counter = 0;
	auto nv = mesh.n_vertices();
	while (ChooseCollapasedEdge(threshold))
	{
		counter++;
		if (nv - counter <= tar_num_v)
		{
			mesh.garbage_collection(true);
			break;
		}
	}
	std::cout << "[V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]" << std::endl;
	std::cout << "Simplify is over!  the iter number is:  " << counter << std::endl;
}

int QEMSimplification::OriginalIndex(int currentid) const
{
	return mesh.property(vindex, mesh.vertex_handle(currentid));
}

const Mesh & QEMSimplification::GetMesh(void) const
{
	return mesh;
}

bool QEMSimplification::ChooseCollapasedEdge(double err_threshhold)
{
	EdgeQEM currentEdge;
	bool done = false;
	int valence = 0;

	do
	{
		done = false;
		if (qpq.size() == 0)
		{
			done = true;
			mesh.garbage_collection(true);
			return false;
		}
		else
		{
			currentEdge = qpq.top();
			qpq.pop();
			valence = mesh.valence(mesh.from_vertex_handle(mesh.halfedge_handle(currentEdge.idx, 0)));
			valence += mesh.valence(mesh.to_vertex_handle(mesh.halfedge_handle(currentEdge.idx, 0)));
		}
	} while ((mesh.status(mesh.edge_handle(mesh.halfedge_handle(currentEdge.idx, 0))).deleted()) || (mesh.property(QuaMat_e, currentEdge.idx) != currentEdge.qem) || (valence > 12));
	if (currentEdge.qem < err_threshhold)
	{
		Mesh::HalfedgeHandle idx = mesh.halfedge_handle(currentEdge.idx, 0);
		std::vector<int> re_in(4);
		re_in[0] = mesh.property(vindex, mesh.from_vertex_handle(idx));
		re_in[1] = mesh.property(vindex, mesh.to_vertex_handle(idx));
		re_in[2] = mesh.property(vindex, mesh.to_vertex_handle(mesh.next_halfedge_handle(idx)));
		re_in[3] = mesh.property(vindex, mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(idx))));
		index_collapse_v_to_v.push_back(re_in);
		OpenMesh::Vec3d& p = mesh.property(new_point, currentEdge.idx);
		mesh.set_point(mesh.to_vertex_handle(idx), p);
		int update_p = mesh.to_vertex_handle(idx).idx();
		mesh.collapse(idx);
		UpdateEdgeCost(update_p);
		return true;
	}
	else
	{
		mesh.garbage_collection(true);
		std::cout << "can not collapsed anymore !   " << "Current minimize qem is : " << currentEdge.qem << std::endl;
		return false;
	}
}

void QEMSimplification::InitialEdgeCost(void)
{
	for (const auto & fh : mesh.faces())
	{
		auto heh = mesh.halfedge_handle(fh);
		const auto & p0 = mesh.point(mesh.from_vertex_handle(heh));
		const auto & p1 = mesh.point(mesh.to_vertex_handle(heh));
		const auto & p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
		auto n = (p1 - p0) % (p2 - p0);
		double area = n.norm();
		mesh.property(area_f, fh) = area;
		n = (area == 0.0) ? Mesh::Point(0) : (n / area);

		double a = n[0];
		double b = n[1];
		double c = n[2];
		double d = -(n | p0);
		Eigen::Matrix4d mat;
		mat << a * a, a * b, a * c, a * d,
			b * a, b * b, b * c, b * d,
			c * a, c * b, c * c, c * d,
			d * a, d * b, d * c, d * d;
		mesh.property(QuaMat_f, fh) = mat;
	}

	for (const auto & vh : mesh.vertices())
	{
		Eigen::Matrix4d mat = Eigen::Matrix4d::Zero();
		for (const auto & vfh : mesh.vf_range(vh))
		{
			mat += mesh.property(QuaMat_f, vfh);
		}
		mesh.property(QuaMat_v, vh) = mat;
	}
	std::vector<EdgeQEM> initialcost;
	initialcost.reserve(mesh.n_edges());
	for (const auto & eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.to_vertex_handle(heh);
		auto vh1 = mesh.from_vertex_handle(heh);
		auto mat = mesh.property(QuaMat_v, vh0) + mesh.property(QuaMat_v, vh1);

		Eigen::Matrix4d mid_m = mat;
		mid_m(3, 0) = mid_m(3, 1) = mid_m(3, 2) = 0.0;
		mid_m(3, 3) = 1.0;
		Eigen::Vector4d b(0.0, 0.0, 0.0, 1.0);
		Eigen::Vector4d v;
		OpenMesh::Vec3d new_p;
		if (fabs(mid_m.determinant()) < 1e-4)
		{
			new_p = (mesh.point(vh0) + mesh.point(vh1)) / 2;
		}
		else
		{
			v = mid_m.colPivHouseholderQr().solve(b);
			new_p[0] = v(0);
			new_p[1] = v(1);
			new_p[2] = v(2);
		}
		mesh.property(new_point, eh) = new_p;
		double err = v.transpose() * mat * v;
		err *= (mesh.valence(vh0) * mesh.valence(vh1));
		err = mesh.is_collapse_ok(heh) ? err : DBL_MAX;
		mesh.property(QuaMat_e, eh) = err;
		initialcost.push_back({ eh, err });
	}
	std::priority_queue<EdgeQEM> qpq_(initialcost.begin(), initialcost.end());
	qpq = qpq_;
}

void QEMSimplification::UpdateEdgeCost(int i)
{
	EdgeQEM e;
	Mesh::VertexHandle v_it = mesh.vertex_handle(i);
	Eigen::Matrix4d q_mat = Eigen::Matrix4d::Zero();
	for (auto vf_it = mesh.vf_begin(v_it); vf_it != mesh.vf_end(v_it); ++vf_it)
	{
		Mesh::FaceVertexIter fv_it = mesh.fv_iter(*vf_it);
		OpenMesh::Vec3d& p = mesh.point(*fv_it);
		++fv_it; OpenMesh::Vec3d& p1 = mesh.point(*fv_it);
		++fv_it; OpenMesh::Vec3d& p2 = mesh.point(*fv_it);
		double area_ = ((p1 - p) % (p2 - p)).norm();
		mesh.property(area_f, *vf_it) = area_;

		OpenMesh::Vec3d& n = ((p1 - p) % (p2 - p)).normalize();
		double a = n[0]; double b = n[1]; double c = n[2]; double d = -(n | p);
		Eigen::Matrix4d mat;
		mat << a * a, a*b, a*c, a*d,
			b*a, b*b, b*c, b*d,
			c*a, c*b, c*c, c*d,
			d*a, d*b, d*c, d*d;
		mesh.property(QuaMat_f, *vf_it) = mat;
		q_mat += mat;
	}
	mesh.property(QuaMat_v, v_it) = q_mat;

	for (auto vv_it = mesh.vv_begin(v_it); vv_it != mesh.vv_end(v_it); ++vv_it)
	{
		q_mat = Eigen::Matrix4d::Zero();
		for (auto vvf_it = mesh.vf_begin(*vv_it); vvf_it != mesh.vf_end(*vv_it); ++vvf_it)
		{
			//q_mat += m.property(area_f,vvf_it)*m.property(QuaMat_f, vvf_it);
			q_mat += mesh.property(QuaMat_f, *vvf_it);
		}
		mesh.property(QuaMat_v, *vv_it) = q_mat;
	}


	for (auto vv_it = mesh.vv_begin(v_it); vv_it != mesh.vv_end(v_it); ++vv_it)
	{
		for (auto e_it = mesh.ve_begin(*vv_it); e_it != mesh.ve_end(*vv_it); ++e_it)
		{
			if (!mesh.property(is_cal_e, *e_it))
			{
				Mesh::HalfedgeHandle he = mesh.halfedge_handle(*e_it, 0);
				Mesh::VertexHandle  v0 = mesh.to_vertex_handle(he);
				Mesh::VertexHandle  v1 = mesh.from_vertex_handle(he);
				Eigen::Matrix4d mat = mesh.property(QuaMat_v, v0) + mesh.property(QuaMat_v, v1);
				Eigen::Matrix4d mid_m = mat;  mid_m(3, 0) = mid_m(3, 1) = mid_m(3, 2) = 0.0; mid_m(3, 3) = 1.0;
				Eigen::Vector4d b(0.0, 0.0, 0.0, 1.0);
				//Eigen::Vector4d v = mid_m.colPivHouseholderQr().solve(b);
				//OpenMesh::Vec3d new_p(v(0), v(1), v(2));
				Eigen::Vector4d v;
				OpenMesh::Vec3d new_p;
				if (fabs(mid_m.determinant()) < 1e-4)
				{
					//cout << mid_m.determinant()<< endl;
					new_p = (mesh.point(v0) + mesh.point(v1)) / 2;
				}
				else
				{
					//v = mid_m.lu().solve(b);
					v = mid_m.colPivHouseholderQr().solve(b);
					new_p[0] = v(0); new_p[1] = v(1); new_p[2] = v(2);
				}
				mesh.property(new_point, *e_it) = new_p;
				double err = v.transpose()*mat*v;
				err *= (mesh.valence(v0) * mesh.valence(v1));
				if (mesh.is_collapse_ok(he) && mesh.is_collapse_ok(mesh.opposite_halfedge_handle(he)))
				{

					mesh.property(QuaMat_e, *e_it) = err;
					e.qem = err;
					e.idx = *e_it;
					qpq.push(e);
					//cout <<m.from_vertex_handle(he).idx() <<"  "<<m.to_vertex_handle(he).idx()<< endl;
				}
				else
				{
					mesh.property(QuaMat_e, *e_it) = DBL_MAX;
				}
				//m.property(edge_l, e_it) = m.calc_edge_length(e_it);
				mesh.property(is_cal_e, *e_it) = true;
			}
		}
	}
	for (auto vv_it = mesh.vv_begin(v_it); vv_it != mesh.vv_end(v_it); ++vv_it)
	{
		for (auto e_it = mesh.ve_begin(*vv_it); e_it != mesh.ve_end(*vv_it); ++e_it)
		{
			mesh.property(is_cal_e, *e_it) = false;
		}
	}
}

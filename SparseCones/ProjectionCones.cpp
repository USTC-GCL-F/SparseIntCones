#include "ProjectionCones.h"


ProjectionCones::ProjectionCones()
{
}

ProjectionCones::~ProjectionCones()
{
}

void ProjectionCones::projection_onto_ori_mesh(const Mesh& mesh, const Mesh& mesh_ori, vector<Opt::VInfo> vinfo_projection)
{
	// ȫ�ֱ���
	// rho_vec; vinfo( û�кϲ� ����ֱ���� )
	// 1.ͶӰ mesh_ori
	//		KdTree
	// 2.������ϲ� mesh_ori
	//		search merge

	std::vector<GG::Triangle> triangles;
	vector<GG::Vec3d> vv;
	vector<unsigned int> idx;
	for (int i = 0; i < mesh_ori.n_vertices(); i++)
	{
		auto v = mesh_ori.vertex_handle(i);
		GG::Vec3d kv;
		auto& p = mesh_ori.point(v);
		kv.x() = p[0];
		kv.y() = p[1];
		kv.z() = p[2];
		vv.push_back(kv);
		idx.push_back(i);
	}
	KdTree kd_tree(vv, idx);

	// �ظ�������Ĳ���
	//sinfo_vec.clear();
	//std::ofstream _out2("projection_before.txt");
	//std::ofstream _out3("projection_after.txt");
	after_projection.clear();
	for (int i = 0; i < vinfo_projection.size(); i++)
	{
		auto v = mesh.vertex_handle(vinfo_projection[i].idx_i);
		auto a_simplify = mesh.point(v);
		auto a = kd_tree.search_nearest_point({ a_simplify[0],a_simplify[1],a_simplify[2] });
		after_projection.push_back(a.id);
		//cout << "ͶӰǰ�� " << vinfo_projection[i].idx_i << endl;
		//_out2 << vinfo_projection[i].idx_i << endl;
		//vinfo_projection[i].idx_i = a.id;
		//cout << "ͶӰ�� " << vinfo_projection[i].idx_i << endl;
		//_out3 << vinfo_projection[i].idx_i << endl;
		
		//SimplifyInfo sinfo;
		//sinfo.idx_i = vinfo_projection[i].idx_i;
		//sinfo.after_reweighed_rho = vinfo_projection[i].after_reweighed_rho;
		//sinfo.gap = vinfo_projection[i].gap;
		//sinfo.after_round= vinfo_projection[i].after_round;
		//auto v = mesh.vertex_handle(vinfo_projection[i].idx_i);
		//auto a = mesh.point(v);
		//sinfo.pi = { a[0],a[1],a[2] };
		//sinfo_vec.push_back(sinfo);
	}
	//_out2.close();
	//_out3.close();
	/*
	for (int i = 0; i < sinfo_vec.size(); i++)
	{
		
		//cones_index.push_back(a.id);
		//cout << "�����񶥵�����: " << v_info_vector[i].idx_i << "   �����ͶӰ --- ԭ���񶥵�����: " << a.id << endl;
	}*/
}

#if 0
void ProjectionCones::projection_onto_ori_mesh_kai(const Mesh& mesh, const Mesh& mesh_ori, vector<KaiCones::VInfo> vinfo_projection)
{
	// ȫ�ֱ���
	// rho_vec; vinfo( û�кϲ� ����ֱ���� )
	// 1.ͶӰ mesh_ori
	//		KdTree
	// 2.������ϲ� mesh_ori
	//		search merge

	std::vector<Triangle> triangles;
	vector<GG::Vec3d> vv;
	vector<unsigned int> idx;
	for (int i = 0; i < mesh_ori.n_vertices(); i++)
	{
		auto v = mesh_ori.vertex_handle(i);
		GG::Vec3d kv;
		auto& p = mesh_ori.point(v);
		kv.x = p[0];
		kv.y = p[1];
		kv.z = p[2];
		vv.push_back(kv);
		idx.push_back(i);
	}
	GG::KdTree kd_tree(vv, idx);

	// �ظ�������Ĳ���
	//sinfo_vec.clear();
	//std::ofstream _out2("projection_before.txt");
	//std::ofstream _out3("projection_after.txt");
	after_projection.clear();
	for (int i = 0; i < vinfo_projection.size(); i++)
	{
		auto v = mesh.vertex_handle(vinfo_projection[i].idx_i);
		auto a_simplify = mesh.point(v);
		auto a = kd_tree.search_nearest_point({ a_simplify[0],a_simplify[1],a_simplify[2] });
		after_projection.push_back(a.id);
		//cout << "ͶӰǰ�� " << vinfo_projection[i].idx_i << endl;
		//_out2 << vinfo_projection[i].idx_i << endl;
		//vinfo_projection[i].idx_i = a.id;
		//cout << "ͶӰ�� " << vinfo_projection[i].idx_i << endl;
		//_out3 << vinfo_projection[i].idx_i << endl;

		//SimplifyInfo sinfo;
		//sinfo.idx_i = vinfo_projection[i].idx_i;
		//sinfo.after_reweighed_rho = vinfo_projection[i].after_reweighed_rho;
		//sinfo.gap = vinfo_projection[i].gap;
		//sinfo.after_round= vinfo_projection[i].after_round;
		//auto v = mesh.vertex_handle(vinfo_projection[i].idx_i);
		//auto a = mesh.point(v);
		//sinfo.pi = { a[0],a[1],a[2] };
		//sinfo_vec.push_back(sinfo);
	}
	//_out2.close();
	//_out3.close();
	/*
	for (int i = 0; i < sinfo_vec.size(); i++)
	{

		//cones_index.push_back(a.id);
		//cout << "�����񶥵�����: " << v_info_vector[i].idx_i << "   �����ͶӰ --- ԭ���񶥵�����: " << a.id << endl;
	}*/
}
#endif
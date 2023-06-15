#include "MeshCut.h"

MeshCut::MeshCut(const Mesh & mesh_)
	:mesh(mesh_)
{
}

MeshCut::~MeshCut(void)
{
}

void MeshCut::mesh_cut(const std::vector<int>& boundary_edge)
{
	if (boundary_edge.empty())
	{
		cuted_mesh = mesh;
		idx_to_mesh_idx.resize(cuted_mesh.n_vertices());
		for (int i = 0; i < mesh.n_vertices(); i++) idx_to_mesh_idx[i] = i;
		return;
	}

	MarkVertexEdge(boundary_edge);
	ComputeValence();

	he_to_idx.clear();
	idx_to_mesh_idx.clear();
	he_to_idx.resize(mesh.n_halfedges());
	Mesh::HalfedgeHandle h_begin;
	for (int i = 0; i < candidate_seam_edge.size(); i++)
	{
		if (candidate_seam_edge[i])
		{
			h_begin = mesh.halfedge_handle(mesh.edge_handle(i), 0);
			break;
		}
	}

	auto h_iter = h_begin;
	int uv_idx = 0;
	do
	{
		he_to_idx[h_iter.idx()] = uv_idx;
		idx_to_mesh_idx.push_back(mesh.to_vertex_handle(h_iter).idx());
		h_iter = mesh.next_halfedge_handle(h_iter);
		while (!candidate_seam_edge[mesh.edge_handle(h_iter).idx()])
		{
			h_iter = mesh.opposite_halfedge_handle(h_iter);
			he_to_idx[h_iter.idx()] = uv_idx;
			if (mesh.is_boundary(h_iter) && !candidate_seam_edge[mesh.edge_handle(mesh.next_halfedge_handle(h_iter)).idx()])
			{
				idx_to_mesh_idx.push_back(mesh.to_vertex_handle(h_iter).idx());
				uv_idx++;
			}
			h_iter = mesh.next_halfedge_handle(h_iter);
		}
		uv_idx++;
	} while (h_iter != h_begin);
	boundary_number = uv_idx;

	for (const auto& vh : mesh.vertices())
	{
		if (candidate_cut_valence[vh.idx()] > 0) continue;
		for (const auto& viheh : mesh.vih_range(vh))
		{
			he_to_idx[viheh.idx()] = uv_idx;
		}
		idx_to_mesh_idx.push_back(vh.idx());
		uv_idx++;
	}

	cuted_mesh.clear();
	for (int i = 0; i < uv_idx; i++)
	{
		cuted_mesh.add_vertex(mesh.point(mesh.vertex_handle(idx_to_mesh_idx[i])));
	}
	for (const auto& fh : mesh.faces())
	{
		std::vector<Mesh::VertexHandle> face_vhandles;
		for (const auto& fheh : mesh.fh_range(fh))
		{
			face_vhandles.push_back(cuted_mesh.vertex_handle(he_to_idx[fheh.idx()]));
		}
		cuted_mesh.add_face(face_vhandles);
	}
}

void MeshCut::mesh_cut(const std::vector<std::vector<int>> & boundary_point)
{
	if (boundary_point.empty())
	{
		cuted_mesh = mesh;
		idx_to_mesh_idx.resize(cuted_mesh.n_vertices());
		for (int i = 0; i < mesh.n_vertices(); i++) idx_to_mesh_idx[i] = i;
		return;
	}

	MarkVertexEdge(boundary_point);
	ComputeValence();
	
	he_to_idx.clear();
	idx_to_mesh_idx.clear();
	he_to_idx.resize(mesh.n_halfedges());
	Mesh::HalfedgeHandle h_begin;
	for (int i = 0; i < candidate_seam_edge.size(); i++)
	{
		if (candidate_seam_edge[i])
		{
			h_begin = mesh.halfedge_handle(mesh.edge_handle(i), 0);
			break;
		}
	}

	auto h_iter = h_begin;
	int uv_idx = 0;
	do
	{
		he_to_idx[h_iter.idx()] = uv_idx;
		idx_to_mesh_idx.push_back(mesh.to_vertex_handle(h_iter).idx());
		h_iter = mesh.next_halfedge_handle(h_iter);
		while (!candidate_seam_edge[mesh.edge_handle(h_iter).idx()])
		{
			h_iter = mesh.opposite_halfedge_handle(h_iter);
			he_to_idx[h_iter.idx()] = uv_idx;
			if (mesh.is_boundary(h_iter) && !candidate_seam_edge[mesh.edge_handle(mesh.next_halfedge_handle(h_iter)).idx()])
			{
				idx_to_mesh_idx.push_back(mesh.to_vertex_handle(h_iter).idx());
				uv_idx++;
			}
			h_iter = mesh.next_halfedge_handle(h_iter);
		}
		uv_idx++;
	} while (h_iter != h_begin);
	boundary_number = uv_idx;

	for (const auto & vh : mesh.vertices())
	{
		if (candidate_cut_valence[vh.idx()] > 0) continue;
		for (const auto & viheh : mesh.vih_range(vh))
		{
			he_to_idx[viheh.idx()] = uv_idx;
		}
		idx_to_mesh_idx.push_back(vh.idx());
		uv_idx++;
	}

	cuted_mesh.clear();
	for (int i = 0; i < uv_idx; i++)
	{
		cuted_mesh.add_vertex(mesh.point(mesh.vertex_handle(idx_to_mesh_idx[i])));
	}
	for (const auto & fh : mesh.faces())
	{
		std::vector<Mesh::VertexHandle> face_vhandles;
		for (const auto & fheh : mesh.fh_range(fh))
		{
			face_vhandles.push_back(cuted_mesh.vertex_handle(he_to_idx[fheh.idx()]));
		}
		cuted_mesh.add_face(face_vhandles);
	}
}

void MeshCut::MarkVertexEdge(const std::vector<int>& boundaryedges)
{
	candidate_seam_vertex.clear();
	candidate_seam_vertex.resize(mesh.n_vertices(), false);
	candidate_seam_edge.clear();
	candidate_seam_edge.resize(mesh.n_edges(), false);

	for (auto e : boundaryedges)
	{
		candidate_seam_edge[e] = true;
		auto eh = OpenMesh::EdgeHandle(e);

		auto vha = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0));
		auto vhb = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 1));
		candidate_seam_vertex[vha.idx()] = true;
		candidate_seam_vertex[vhb.idx()] = true;
	}
}

void MeshCut::MarkVertexEdge(const std::vector<std::vector<int>>& boundarypoints)
{
	candidate_seam_vertex.clear();
	candidate_seam_vertex.resize(mesh.n_vertices(), false);
	candidate_seam_edge.clear();
	candidate_seam_edge.resize(mesh.n_edges(), false);

	for (const auto & pts : boundarypoints)
	{
		candidate_seam_vertex[pts[0]] = true;
		for (size_t i = 0; i + 1 < pts.size(); i++)
		{
			candidate_seam_vertex[pts[i + 1]] = true;
			for (const auto& heh : mesh.voh_range(mesh.vertex_handle(pts[i])))
			{
				if (mesh.to_vertex_handle(heh) == mesh.vertex_handle(pts[i + 1]))
				{
					candidate_seam_edge[heh.idx() >> 1] = true;
				}
			}
		}
	}
}

void MeshCut::ComputeValence(void)
{
	candidate_cut_valence.clear();
	candidate_cut_valence.resize(mesh.n_vertices(), 0);
	for (const auto & eh : mesh.edges())
	{
		if (candidate_seam_edge[eh.idx()])
		{
			auto heh = mesh.halfedge_handle(eh, 0);
			++candidate_cut_valence[mesh.from_vertex_handle(heh).idx()];
			++candidate_cut_valence[mesh.to_vertex_handle(heh).idx()];
		}
	}
}

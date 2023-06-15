#include "MinimalSpanningTree.h"

MinimalSpanningTree::~MinimalSpanningTree(void)
{
}

MinimalSpanningTree::MinimalSpanningTree(const Mesh & m)
	:mesh(m),
	nv(0),
	ne(0),
	findgeodesic(mesh)
{
}

MinimalSpanningTree::MinimalSpanningTree(const Mesh & m, const std::vector<std::vector<double>> & l)
	: MinimalSpanningTree(m)
{
	findgeodesic.SetEdgeWeight(l);
}

MinimalSpanningTree::MinimalSpanningTree(const Mesh & m, const std::vector<int> & vertex)
	:MinimalSpanningTree(m, vertex, {})
{
}

MinimalSpanningTree::MinimalSpanningTree(const Mesh & m, const std::vector<int> & vertex, const std::vector<int> & isaccessible)
	:mesh(m),
	nodes(vertex),
	nv(static_cast<int>(vertex.size())),
	ne(nv * (nv - 1) / 2),
	findgeodesic(mesh)
{
	shortestedges.reserve(ne);
	weight.reserve(ne);
	std::vector<int> is_target(mesh.n_vertices(), 0);
	for (int i = 1; i < nv; i++)
	{
		is_target[nodes[i - 1]] = 1;
		findgeodesic.Compute(nodes[i], is_target, isaccessible, i);
		for (int j = 0; j < i; j++)
		{
			shortestedges.push_back(findgeodesic.GetPath(nodes[j]));
			weight.push_back(findgeodesic.GetLength(nodes[j]));
		}
	}
}

void MinimalSpanningTree::SetNodes(const std::vector<int>& vertex, const std::vector<int> & isaccessible)
{
	nodes = vertex;
	nv = static_cast<int>(vertex.size());
	ne = nv * (nv - 1) / 2;
	shortestedges.clear();
	shortestedges.reserve(ne);
	weight.clear();
	weight.reserve(ne);
	std::vector<int> is_target(mesh.n_vertices(), 0);
	for (int i = 1; i < nv; i++)
	{
		is_target[nodes[i - 1]] = 1;
		findgeodesic.Compute(nodes[i], is_target, isaccessible, i);
		for (int j = 0; j < i; j++)
		{
			shortestedges.push_back(findgeodesic.GetPath(nodes[j]));
			weight.push_back(findgeodesic.GetLength(nodes[j]));
		}
	}
}

void MinimalSpanningTree::Kruskal(void)
{
	spanningtreeedges.clear();
	std::vector<int> spanning_tree_current(nv);
	for (int i = 0; i < nv; ++i)
	{
		spanning_tree_current[i] = i;
	}
	std::vector<std::pair<int, int>> edges = GetEdges(); // edges in the graph
	int index = 1;
	int j = 0;
	while (index < nv)
	{
		auto m = spanning_tree_current[edges[j].first];
		auto n = spanning_tree_current[edges[j].second];
		if (m != n)
		{
			spanningtreeedges.push_back(shortestedges[edges[j].second + edges[j].first * (edges[j].first - 1) / 2]);
			index++;
			for (int i = 0; i < nv; ++i)
			{
				if (spanning_tree_current[i] == n)
				{
					spanning_tree_current[i] = m;
				}
			}
		}
		j++;
	}
}

void MinimalSpanningTree::Prim(void)
{
	spanningtreeedges.clear();
	std::vector<std::pair<int, int>> edges = GetEdges();// edges in the graph
	spanningtreeedges.push_back(shortestedges[edges[0].second + edges[0].first * (edges[0].first - 1) / 2]);
	std::vector<int> istarget(mesh.n_vertices(), 0);
	for (const auto & s : spanningtreeedges[0])
	{
		istarget[s] = 1;
	}

	std::vector<int> source;
	for (const auto & i : nodes)
	{
		if (i != nodes[edges[0].first] && i != nodes[edges[0].second])
		{
			source.push_back(i);
		}
	}
	while (!source.empty())
	{
		double minlen = DBL_MAX;
		std::vector<int> minpath;
		int minsource;
		for (const auto & souv : source)
		{
			auto tarv = findgeodesic.Compute(souv, istarget, {});
			auto len = findgeodesic.GetLength(tarv);
			if (len < minlen)
			{
				minlen = len;
				minpath = findgeodesic.GetPath(tarv);
				minsource = souv;
			}
		}
		spanningtreeedges.push_back(minpath);
		for (const auto & vid : minpath)
		{
			istarget[vid] = 1;
		}
		source.erase(std::find(source.begin(), source.end(), minsource));
	}
}

void MinimalSpanningTree::PrimAlter(void)
{
	spanningtreeedges.clear();
	std::vector<std::pair<int, int>> edges = GetEdges(); // edges in the graph
	spanningtreeedges.push_back(shortestedges[edges[0].second + edges[0].first * (edges[0].first - 1) / 2]);
	std::vector<int> istarget(mesh.n_vertices(), 0);
	std::vector<int> isaccessible(mesh.n_vertices(), 1);
	for (const auto & s : spanningtreeedges[0])
	{
		isaccessible[s] = 0;
	}
	isaccessible[nodes[edges[0].first]] = 1;
	isaccessible[nodes[edges[0].second]] = 1;
	istarget[nodes[edges[0].first]] = 1;
	istarget[nodes[edges[0].second]] = 1;

	std::vector<int> source;
	for (const auto & i : nodes)
	{
		if (i != nodes[edges[0].first] && i != nodes[edges[0].second])
		{
			source.push_back(i);
		}
	}
	while (!source.empty())
	{
		double minlen = DBL_MAX;
		std::vector<int> minpath;
		int minsource;
		for (const auto & souv : source)
		{
			auto tarv = findgeodesic.Compute(souv, istarget, isaccessible);
			auto path = findgeodesic.GetPath(tarv);
			auto len = findgeodesic.GetLength(tarv);
			if (len < minlen)
			{
				minlen = len;
				minpath = path;
				minsource = souv;
			}
		}
		spanningtreeedges.push_back(minpath);
		for (const auto & vid : minpath)
		{
			isaccessible[vid] = 0;
		}
		isaccessible[minpath.front()] = 1;
		isaccessible[minpath.back()] = 1;
		istarget[minpath.front()] = 1;
		istarget[minpath.back()] = 1;

		std::vector<int> newsource;
		for (const auto & souv : source)
		{
			if (souv != minsource)
			{
				newsource.push_back(souv);
			}
		}
		source = newsource;
	}
}

const std::vector<std::vector<int>> & MinimalSpanningTree::GetSpanningTreeEdges(void) const
{
	return spanningtreeedges;
}

std::vector<std::pair<int, int>> MinimalSpanningTree::GetEdges(void) const
{
	std::vector<std::pair<int, int>> edges;
	edges.reserve(ne);
	int k = 0;
	for (int i = 1; i < nv; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			edges.push_back({ i, j });
			k++;
		}
	}
	std::sort(edges.begin(), edges.end(), [&](const std::pair<int, int> & lhs, const std::pair<int, int> & rhs) { return weight[lhs.second + lhs.first * (lhs.first - 1) / 2] < weight[rhs.second + rhs.first * (rhs.first - 1) / 2]; });
	return edges;
}

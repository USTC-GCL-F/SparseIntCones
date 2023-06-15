#include "ShortestPath.h"
#include <queue>

ShortestPath::ShortestPath(const Mesh & mesh_)
	:mesh(mesh_),
	dist(mesh_.n_vertices()),
	prev(mesh_.n_vertices()),
	G(mesh_.n_vertices())
{
	for (const auto& vh : mesh_.vertices())
	{
		auto & vecpair = G[vh.idx()];
		const auto& p0 = mesh_.point(vh);
		for (const auto& vvh : mesh_.vv_range(vh))
		{
			vecpair.push_back({ vvh.idx(), (mesh_.point(vvh) - p0).norm() });
		}
	}
}

ShortestPath::~ShortestPath(void)
{
}

void ShortestPath::SetEdgeWeight(const std::vector<std::vector<double>> & weight)
{
	G.clear();
	G.resize(mesh.n_vertices());
	for (const auto & vh : mesh.vertices())
	{
		auto vid = vh.idx();
		auto & vecpair = G[vid];
		int counter = 0;
		for (const auto & vvh : mesh.vv_range(vh))
		{
			vecpair.push_back({ vvh.idx(), weight[vid][counter] });
			counter++;
		}
	}
}

int ShortestPath::Compute(const int & source, const std::vector<int>& istarget, const std::vector<int> & isaccessible, const int& size)
{
	int counter = 0;
	int closest = -1;
	auto hasaccessible = !isaccessible.empty();
	dist.assign(mesh.n_vertices(), DBL_MAX);
	prev.assign(mesh.n_vertices(), -1);
	dist[source] = 0;
	auto cmp = [](const std::pair<int, double> & lhs, const std::pair<int, double> & rhs) { return lhs.second > rhs.second; };
	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, decltype(cmp)> Q(cmp);
	Q.push(std::pair<int, double>(source, 0.0));
	while (!Q.empty())
	{
		const auto & top = Q.top();
		int i = top.first;
		double c = top.second;
		Q.pop();
		
		if (c <= dist[i])
		{
			if (istarget[i])
			{
				if (counter == 0)
				{
					closest = i;
				}
				counter++;
				if (counter == size)
				{
					break;
				}
			}
			for (const auto & pi : G[i])
			{
				int i2 = pi.first;
				double c2 = pi.second;
				if ((c2 + dist[i] < dist[i2]) && (!hasaccessible || isaccessible[i2]))
				{
					dist[i2] = c2 + dist[i];
					prev[i2] = i;
					Q.push({ i2, dist[i2] });
				}
			}
		}
	}
	return closest;
}

std::vector<int> ShortestPath::GetPath(const int & tar_v) const
{
	std::vector<int> geo;
	int tmp = tar_v;
	do
	{
		geo.push_back(tmp);
		tmp = prev[tmp];
	} while (tmp != -1);
	return geo;
}

const double & ShortestPath::GetLength(const int & target) const
{
	return dist[target];
}

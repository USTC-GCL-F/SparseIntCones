#pragma once
#include <string>
#include <fstream>
#include "MeshCache.h"

class Algorithm
{
	Algorithm();
	~Algorithm();
public:
	static void Dijkstra_all(MeshCache& MCache, int k);
	static void Dijkstra_group(MeshCache& MCache, std::vector<int>& lmk);
	static void Kruskal(MeshCache& MCache, std::vector<int>& lmk, std::vector<int>& cutvertex, std::vector<int>& cutedge);
	static void Kruskal(std::vector<int>& lmk, std::priority_queue<PathInfo> que, std::vector<PathInfo>& Result);
	static void FindPath(std::vector<int>& v_p, int e_p, std::vector<int>& path_v);
	static void UpdateNeighbourhood(MeshCache& MCache, int k, int v);
	static void Dijkstra_with_restrict(MeshCache& MCache, int s_p, std::vector<double>& weight, std::vector<int>&, std::vector<double>&);
	static void Dijkstra_with_nearest2(MeshCache& MCache, int s_p, std::vector<int>& is_target, std::vector<int>& path);
};



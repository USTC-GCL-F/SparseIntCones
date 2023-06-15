#pragma once
#include "Mesh/MeshDefinition.h"
#include "ShortestPath.h"

// Class of computing minimal spanning tree on a mesh
class MinimalSpanningTree
{
public:
	~MinimalSpanningTree(void);
	MinimalSpanningTree(const Mesh & m);
	MinimalSpanningTree(const Mesh & m, const std::vector<std::vector<double>> & l);
	MinimalSpanningTree(const Mesh & m, const std::vector<int> & vertex);
	MinimalSpanningTree(const Mesh & m, const std::vector<int> & vertex,
		const std::vector<int> & isaccessible);

	// Set the nodes of the graph
	void SetNodes(const std::vector<int> & vertex, const std::vector<int> & isaccessible = {});

	// Kruskal's algorithm to compute minimal spanning tree
	// Allow duplicate paths
	void Kruskal(void);

	// Prim's algorithm
	// Compute a new edge connected to any vertex on the existing tree
	void Prim(void);

	// Prim's algorithm
	// Compute a new edge connected to nodes and prevent duplicate path
	void PrimAlter(void);

	// Get the final result
	const std::vector<std::vector<int>> & GetSpanningTreeEdges(void) const;

private:
	// Get sorted edges (i, j) (1 <= j < i <= nv - 1)
	std::vector<std::pair<int, int>> GetEdges(void) const;
	
	const Mesh & mesh; // constant reference to the mesh
	std::vector<int> nodes;  // indices of source nodes
	int nv;                  // number of nodes
	int ne;                  // number of edges
	std::vector<double> weight; // adjacent length
	std::vector<std::vector<int>> shortestedges; // adjacent path
	ShortestPath findgeodesic; // shortest path calculator
	std::vector<std::vector<int>> spanningtreeedges; // result
};

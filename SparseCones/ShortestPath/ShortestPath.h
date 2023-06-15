#pragma once
#include "Mesh/MeshDefinition.h"

// Class of shortest path computation on a mesh
class ShortestPath
{
public:
	// Construction from a mesh
	ShortestPath(const Mesh & mesh_);

	// Destruction
	virtual ~ShortestPath(void);

	// Set custom non-negative weight of edges
	// Caution: weight must have same structure with the mesh one-ring
	void SetEdgeWeight(const std::vector<std::vector<double>> & weight);
	
	// Compute the shortest path from a source vertex to a target vertex set
	// Return the index of the closest target vertex
	// source: index of the source vertex
	// istarget: if a vertex is set to be target (1) or not (0)
	// isaccessible: if a vertex is accessible (1) or not (0)
	// size: number of target need to be computed
	int Compute(const int& source, const std::vector<int>& istarget,
		const std::vector<int>& isaccessible, const int& size = 1);
	
	// Get the shortest path to a target vertex
	// This function has to be called after the shortest path is computed
	std::vector<int> GetPath(const int & target) const;

	// Get the length of the shortest path to a target vertex (reversed)
	// This function has to be called after the shortest path is computed
	const double & GetLength(const int & target) const;

private:
	const Mesh & mesh; // constant reference to a mesh
	std::vector<double> dist; // distance from a source to each vertex
	std::vector<int> prev; // previous vertex of each vertex
	std::vector<std::vector<std::pair<int, double>>> G;// adjacent matrix
};

#pragma once
#include "Mesh/MeshDefinition.h"

struct EdgeInfo
{
	//	int edge_id;
	Mesh::HalfedgeHandle idx;
	double ecost;
	inline    friend bool      operator<(const EdgeInfo & lhs, const EdgeInfo & rhs) { return (lhs.ecost < rhs.ecost); }
	inline    friend bool      operator>(const EdgeInfo & lhs, const EdgeInfo & rhs) { return (lhs.ecost > rhs.ecost); }
};

struct EdgeInfo_
{
	//	int edge_id;
	Mesh::HalfedgeHandle idx;
	double ecost;
	inline    friend bool      operator<(const EdgeInfo_ & lhs, const EdgeInfo_ & rhs) { return (lhs.ecost > rhs.ecost); }
	inline    friend bool      operator>(const EdgeInfo_ & lhs, const EdgeInfo_ & rhs) { return (lhs.ecost < rhs.ecost); }
};

struct EdgeQEM
{
	Mesh::EdgeHandle idx;
	double qem;
	inline    friend bool      operator<(const EdgeQEM & lhs, const EdgeQEM & rhs) { return (lhs.qem > rhs.qem); }
	inline    friend bool      operator>(const EdgeQEM & lhs, const EdgeQEM & rhs) { return (lhs.qem < rhs.qem); }
};

class General_Struct
{
public:
	General_Struct();
	~General_Struct();
};


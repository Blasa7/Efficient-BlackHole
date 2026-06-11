#include <Util.hpp>
#include <set>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include <algorithm>

#pragma once

template <uint32_t dimensions>
class node{
protected :
	int nodeId;
	double* points;

public :
	bool setPosition(int nLabel, double* ar);
	double getValue(int idx);
	void setValue(double value,int idx);
	int getID() const;
};

// Blackhole node is a node which manage neighbor information
template <uint32_t dimensions>
class blackHoleNode : public node<dimensions> {
private :
	int clusterId;
	int degree;
	// Set of edges
	std::vector<int> eSet;
public :
	blackHoleNode(int nNodeId, int other);
	~blackHoleNode();
	void setClusterId(int x);
	int getClusterId();
	bool isnot_labeled();
	bool setEdge(int other);
	int getDegree();
	double* getValues();
	void setDegree(int x);
	std::vector<int>* getEdgeSet();
	bool findEdge(int origin, int id);
	inline bool operator<(const blackHoleNode<dimensions>& i) const{
		return this->nodeId < i.getID();	
	};
};

#include "blackHoleNode.cpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <map>
#include <vector>
#include <exponentVar.hpp>
#include <octTree.hpp>

#pragma once

template <uint32_t dimensions>
struct by_id { 
    bool operator()(blackHoleNode<dimensions>* const &a, blackHoleNode<dimensions>* const &b) const {
        return a->getID() < b->getID();
    }
};

template <uint32_t dimensions>
class nodeCollection{
private :
	std::map<int, blackHoleNode<dimensions>*> nodeMap;
	std::vector<blackHoleNode<dimensions>*> nodeVec;
	int* degMat;
	float** adjMat;
public :
	void putNode(int nNodeId, int other);
	void clearAll();
	void copyToVector();
	bool checkEdge(int id1_notMinus, int id2_notMinus);
	void setAdjMat(int maxValue);
	void setDegMat(int maxValue);
	int getSumOfDegree();
	bool checkExists(int nID);
	void degreeSet();
	void clearClusterId();
	std::vector<blackHoleNode<dimensions>*>* getNodeVec();
	std::map<int, blackHoleNode<dimensions>*>& getNodeMap();

	//Barneshut
	double addRepulsionDir(blackHoleNode<dimensions>* unp, double* dir, exponentVar& expVar,OctTree<dimensions>* octTree);
	double addAttractionDirA(blackHoleNode<dimensions>* unp, double* dir, exponentVar& expVar,OctTree<dimensions>* octTree);
	void setDir(blackHoleNode<dimensions>* unp, double* dir,  exponentVar& expVar,OctTree<dimensions>* octTree);
	double getEnergy(blackHoleNode<dimensions>* unp, exponentVar& expVar,OctTree<dimensions>* octTree);
	double getEnergyR(blackHoleNode<dimensions>* unp , exponentVar& expVar,OctTree<dimensions>* octTree);
	double getEnergyAA(blackHoleNode<dimensions>* unp , exponentVar& expVar, OctTree<dimensions>* octTree);
	double findInitEnergy(exponentVar& expVar, OctTree<dimensions>* octTree);
};

#include "nodeCollection.cpp"

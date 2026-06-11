#pragma once

#include "Util.hpp"
#include "blackHoleNode.hpp"


template <uint32_t dimensions>
blackHoleNode<dimensions>::~blackHoleNode(){
	delete this->points;
}

template <uint32_t dimensions>
bool node<dimensions>::setPosition(int nI, double* ar){
	nodeId = nI;
	points = new double[dimensions];
	for (int temp = 0; temp < dimensions; temp++){
		points[temp] = ar[temp];
	}
	return true;
}

template <uint32_t dimensions>
int node<dimensions>::getID() const{
	return nodeId;
}

template <uint32_t dimensions>
double node<dimensions>::getValue(int idx){
	return points[idx];
}

template <uint32_t dimensions>
void node<dimensions>::setValue(double value, int idx){
	points[idx] = value;
}


template <uint32_t dimensions>
double* blackHoleNode<dimensions>::getValues(void){
	return this->points;
}

template <uint32_t dimensions>
void blackHoleNode<dimensions>::setDegree(int x){
	degree = x;
}

template <uint32_t dimensions>
void blackHoleNode<dimensions>::setClusterId(int x){
	clusterId = x;
}

template <uint32_t dimensions>
int blackHoleNode<dimensions>::getClusterId(){
	return clusterId;
}

template <uint32_t dimensions>
bool blackHoleNode<dimensions>::isnot_labeled(){
	if(clusterId == -1){
		return false;
	}
	return true;
}

template <uint32_t dimensions>
std::vector<int>* blackHoleNode<dimensions>::getEdgeSet(){
	return &eSet;
}

template <uint32_t dimensions>
bool blackHoleNode<dimensions>::setEdge(int e){
	eSet.push_back(e);
	return true;
}

template <uint32_t dimensions>
bool blackHoleNode<dimensions>::findEdge(int origin, int id){
	return std::binary_search(eSet.begin(), eSet.end(), id);
}

template <uint32_t dimensions>
blackHoleNode<dimensions>::blackHoleNode(int nNodeId, int ep){
	this->nodeId = nNodeId;
	this->points = new double[dimensions];
	for (int i = 0; i < dimensions; i++){
		this->points[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX)-0.5f;
	}
	clusterId = -1;
	degree = 1;
	eSet.push_back(ep);
}


template <uint32_t dimensions>
int blackHoleNode<dimensions>::getDegree(){
	return degree;
}

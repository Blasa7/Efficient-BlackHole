#pragma once

#include "nodeCollection.hpp"



template <uint32_t dimensions>
void nodeCollection<dimensions>::putNode(int nNodeId, int other){
	if(degMat[nNodeId] == 0){	//not exists!
		degMat[nNodeId]++;
		blackHoleNode<dimensions>* newNode = new blackHoleNode<dimensions>(nNodeId + 1, other + 1);
		nodeMap[nNodeId + 1] = newNode;
	}
	else{	
		if(nodeMap[nNodeId + 1]->findEdge(nNodeId + 1, other + 1) == false){	
			degMat[nNodeId]++;
			nodeMap[nNodeId + 1]->setEdge(other + 1);
		}
	}
}

template <uint32_t dimensions>
void nodeCollection<dimensions>::copyToVector(){

	for(typename std::map<int, blackHoleNode<dimensions>*>::iterator it = nodeMap.begin(); it != nodeMap.end(); it++){
		nodeVec.push_back(it->second);
	}
	std::sort(nodeVec.begin(), nodeVec.end(), by_id<dimensions>());
}

template <uint32_t dimensions>
void nodeCollection<dimensions>::setAdjMat(int maxValue){
	int len = maxValue;
	adjMat = new float*[len];
	for(int i = 0; i < len; i++){
		adjMat[i] = new float[len];
		//memset(adjMat[i], 0, sizeof(float)*len);
		for(int k = 0; k < len; k++){
			adjMat[i][k] = 0;
		}
	}
}

template <uint32_t dimensions>
void nodeCollection<dimensions>::setDegMat(int maxValue){
	int len = maxValue;
	degMat = new int[len];
	for(int i = 0; i < len; i++){
		degMat[i] = 0;
	}
}

template <uint32_t dimensions>
bool nodeCollection<dimensions>::checkEdge(int id1_notMinus, int id2_notMinus){
	return adjMat[id1_notMinus-1][id2_notMinus-1] != 0;
}

template <uint32_t dimensions>
int nodeCollection<dimensions>::getSumOfDegree(){
	std::vector<blackHoleNode<dimensions>*>* vect = &nodeVec;
	int s = 0;
	for(unsigned int i = 0; i < nodeVec.size(); i++){
		s += (*vect)[i]->getDegree();
	}
	return s;
}


template <uint32_t dimensions>
double nodeCollection<dimensions>::findInitEnergy(exponentVar& expVar, OctTree<dimensions>* octTree){
	
	std::vector<blackHoleNode*>* vect = &nodeVec;
	double s = 0.0;

	for(unsigned int i = 0; i < nodeVec.size(); i++){
		s += getEnergy((*vect)[i], expVar, octTree);
	}
	return s;
}

template <uint32_t dimensions>
void nodeCollection<dimensions>::degreeSet(){
	std::vector<blackHoleNode<dimensions>*>* vect = &nodeVec;
	for(unsigned int i = 0; i < nodeVec.size(); i++){
		(*vect)[i]->setDegree(degMat[i]);
	}
}

template <uint32_t dimensions>
void nodeCollection<dimensions>::clearAll(){
	delete [] degMat;

	for(unsigned int i = 0; i < nodeMap.size(); i++){
		delete nodeMap[i];
	}


	nodeMap.clear();
	nodeVec.clear();
}

template <uint32_t dimensions>
std::vector<blackHoleNode<dimensions>*>* nodeCollection<dimensions>::getNodeVec(){
	return &nodeVec;
}

template <uint32_t dimensions>
std::map<int, blackHoleNode<dimensions>*>& nodeCollection<dimensions>::getNodeMap(){
	return nodeMap;
}
/*Use*/
template <uint32_t dimensions>
double nodeCollection<dimensions>::getEnergyR(blackHoleNode<dimensions>* unp, exponentVar& expVar, OctTree<dimensions>* tree){

	if(tree == NULL || unp->getDegree() == 0 )
		return 0.0;
		
	if(tree->node == unp)
		return 0.0;	

	double repuExponent = expVar.getRepuExponent();	
	double repuFactor = expVar.getRepuFactor();

	double* pos = unp->getValues();
	
	double* pos2 = tree->getValues();

	double treeWidth = tree->getWidth();
	double dist = CalculateUtil::calcDist_DIM<dimensions>(pos, pos2);
	if (dist == 0.0) return 0.0;

	if ( tree->childCount > 0 && dist < 1.0 * treeWidth){
		double energy = 0.0;
		for(int i = 0; i < tree->getLength(); i++){
				energy += getEnergyR(unp, expVar, tree->children[i]);
				//std::cout<<energy<<"\t";
		}
		return energy;
	}

	if (dist == 0.0) return 0.0;

	if(repuExponent == 0.0){
		return -repuFactor * unp->getDegree() * tree->weight * log(dist);
	}
	else{
		return -repuFactor * unp->getDegree() * tree->weight * pow(dist, repuExponent) / repuExponent;
	}
}

/*Use*/
template <uint32_t dimensions>
double nodeCollection<dimensions>::addRepulsionDir(blackHoleNode<dimensions>* unp, double* dir, exponentVar& expVar, OctTree<dimensions>* tree){

	//if(tree==NULL || tree->node == NULL ||tree->node->getID() == unp->getID() || unp->getDegree() == 0)
	if(tree==NULL ||  tree->node == unp)
		return 0.0;

	if(unp->getDegree()==0){
		return 0.0;
	}


	//std::cout<<"node id = "<<unp->getID()<<"  Tree node id = "<<tree->node->getID()<<std::endl;

	double repuFactor = expVar.getRepuFactor();
	double repuExponent = expVar.getRepuExponent();


	double* pos = unp->getValues();
	double* pos2 = tree->getValues();
	double dist = CalculateUtil::calcDist_DIM<dimensions>(pos, pos2);


	if (dist == 0.0) return 0.0;

	if (tree->childCount > 0 && dist < 1.0 * tree->getWidth()){
		for(int i = 0; i < tree->childLength; i++){
			addRepulsionDir(unp, dir, expVar, tree->children[i]);
		}
		return 0.0;
	}

	double tmp = repuFactor * unp->getDegree() * tree->getWeight() * pow(dist, repuExponent-2);

	for (int z = 0; z < dimensions; z++){
		dir[z] -= (pos2[z] - pos[z]) * tmp;
	}
	return tmp * abs(repuExponent-1);
}	

/*Use*/
template <uint32_t dimensions>
double nodeCollection<dimensions>::getEnergyAA(blackHoleNode<dimensions>* unp, exponentVar& expVar, OctTree<dimensions>* tree){
	
	double attrExponent = expVar.getAttrExponent();
	std::vector<blackHoleNode<dimensions>*>* vect = &nodeVec;
	
	double* pos = unp->getValues();

	std::vector<int>* edgeList = unp->getEdgeSet();

	double energy = 0.0;

	for(unsigned int i = 0; i < (*edgeList).size(); i++){
		int value = (*edgeList)[i]-1;
		double* pos2 = (*vect)[value]->getValues();
		double dst = CalculateUtil::calcDist_DIM<dimensions>(pos, pos2);
		if (attrExponent == 0.0){
			energy += log(dst);
		}
		else{
			energy += pow(dst, attrExponent) / attrExponent;
		}
	}
	return energy;
}


/*Use*/
template <uint32_t dimensions>
double nodeCollection<dimensions>::addAttractionDirA(blackHoleNode<dimensions>* unp, double* dir, exponentVar& expVar, OctTree<dimensions>* tree){
	
	std::vector<blackHoleNode<dimensions>*>* vect = &nodeVec;
	
	double attrExponent = expVar.getAttrExponent();

	if(unp->getDegree() == 0){
		return 0.0;
	}
	double* pos = unp->getValues();

	std::vector<int>* edgeList = unp->getEdgeSet();
	
	for(unsigned int i = 0; i <  (*edgeList).size(); i++){
		int value = (*edgeList)[i]-1;
		double* pos2 = (*vect)[value]->getValues();
		double dist = CalculateUtil::calcDist_DIM<dimensions>(pos, pos2);
		if(dist == 0.0){
			continue;
		}
		double tmp = pow(dist, attrExponent-2);
		
		for (int z = 0; z < dimensions; z++){
			dir[z] += (pos2[z] - pos[z]) * tmp;
		}
	}
	return 0.0;
}


/*Use*/
template <uint32_t dimensions>
double nodeCollection<dimensions>::getEnergy(blackHoleNode<dimensions>* unp, exponentVar& expVar, OctTree<dimensions>* tree){
	double gr = getEnergyR(unp, expVar, tree);
	double ga = getEnergyAA(unp, expVar, tree);
	
	return gr + ga;
}

/*Use*/
template <uint32_t dimensions>
void nodeCollection<dimensions>::setDir(blackHoleNode<dimensions>* unp, double* dir, exponentVar& expVar, OctTree<dimensions>* tree){
	addRepulsionDir(unp, dir, expVar, tree);
	addAttractionDirA(unp, dir, expVar, tree);
}

template <uint32_t dimensions>
void nodeCollection<dimensions>::clearClusterId(){
	std::vector<blackHoleNode*>* vect = &nodeVec;
	for (unsigned int i = 0; i < nodeVec.size(); i++){
		(*vect)[i]->setClusterId(-1);
	}
}

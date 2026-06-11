#ifndef _OCT_H_
#define _OCT_H_

#include "Util.hpp"
#include "blackHoleNode.hpp"
#include "memoryManager.hpp"

template <uint32_t dimensions>
class memoryManager;

template <uint32_t dimensions>
class OctTree{
public :
	static const int MAX_DEPTH = 18;
	blackHoleNode<dimensions>* node;
	OctTree<dimensions>** children;
	int childCount;
	int childLength;
	double* position;
	double weight;
	double* minPos;
	double* maxPos;

public :

	OctTree();

	~OctTree();

	OctTree(blackHoleNode<dimensions>* node, double* position, double* minPos, double* maxPos);

	void setElement(blackHoleNode<dimensions>* node, double* position, double* minPos, double* maxPos, memoryManager<dimensions>* mgr);

	double getPosX();

	double getPosY();

	double* getValues();

	blackHoleNode<dimensions>* getNode();

	int getLength();

	void LengthIncrease();

	double getWeight();

	bool getUsed();

	void addNode(blackHoleNode<dimensions>* newNode, double* newPos, int depth, memoryManager<dimensions>* mgr);

	void addNode2(blackHoleNode<dimensions>* newNode, double* newPos, int depth, memoryManager<dimensions>* mgr, bool st);

	double getWidth();

	int getHeight();

	void removeNode(blackHoleNode<dimensions>* oldNode, double* oldPos, int depth, memoryManager<dimensions>* mgr);

	void clearMemory();
	 
	void clearMemory(OctTree<dimensions>* p);
};

#include "octTree.cpp"

#endif

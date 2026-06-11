
#ifndef _MEM_H_
#define _MEM_H_

#include "octTree.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>


template <uint32_t dimensions>
class OctTree;

template <uint32_t dimensions>
class memoryManager{
private :
	OctTree<dimensions>* temp;
	std::vector<OctTree<dimensions>**> childVec;	//new OctTree*[8]
	int current;
	int child_current;
	int prior;
	int maxNum;

public :

	memoryManager(int nodeNum);

	OctTree<dimensions>* get_Instance();

	OctTree<dimensions>** get_children();

	void dealloc();

	int getCurrent();

	int getChildCurrent();

	void setCurrent(int x);

	void setChildCurrent(int x);

	void takeAPicture();

	void restore();

	void swap(OctTree<dimensions>* swapper);

};

#include "memoryManager.cpp"

#endif

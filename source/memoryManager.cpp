#pragma once

#include "memoryManager.hpp"

template <uint32_t dimensions>
memoryManager<dimensions>::memoryManager(int nodeNum){
	current = 0;
	child_current = 0;
	maxNum = nodeNum * 12;
	prior = -1;
	temp = new OctTree<dimensions>[maxNum];
	
	for(int j = 0; j < maxNum; j++){
		childVec.push_back(new OctTree<dimensions>*[(int)pow(2.0, dimensions)]);
	}
}

template <uint32_t dimensions>
OctTree<dimensions>** memoryManager<dimensions>::get_children(){
	return childVec[child_current++];
}

template <uint32_t dimensions>
OctTree<dimensions>* memoryManager<dimensions>::get_Instance(){
	return &(temp[current++]);
}

template <uint32_t dimensions>
int memoryManager<dimensions>::getCurrent(){
	return current;
}

template <uint32_t dimensions>
int memoryManager<dimensions>::getChildCurrent(){
	return child_current;
}

template <uint32_t dimensions>
void memoryManager<dimensions>::setCurrent(int x){
	current = x;
}

template <uint32_t dimensions>
void memoryManager<dimensions>::setChildCurrent(int x){
	child_current = x;
}

template <uint32_t dimensions>
void memoryManager<dimensions>::dealloc(){

	delete[] temp;

	for(int i = 0; i < maxNum; i++){
		delete childVec[i];
	}
	std::cout << "childVec removed!" << std::endl;

	childVec.clear();
}

template <uint32_t dimensions>
void memoryManager<dimensions>::takeAPicture(){
	prior = current;
}

template <uint32_t dimensions>
void memoryManager<dimensions>::restore(){
	current = prior;
	child_current = prior;
}

template <uint32_t dimensions>
void memoryManager<dimensions>::swap(OctTree<dimensions>* swapper){
	if(prior == -1 && current == prior)
		return;
	prior++;
	OctTree<dimensions> s = *swapper;
	*swapper = temp[prior];
	temp[prior] = s;

	OctTree<dimensions>** ss = swapper->children;
	swapper->children = temp[prior].children;
	temp[prior].children = ss;
}

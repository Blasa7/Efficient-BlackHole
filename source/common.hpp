#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Some typedefs to help with readability.
typedef std::vector<std::pair<uint32_t, uint32_t>> EdgeList;
typedef std::vector<std::vector<uint32_t>> AdjacencyList;
typedef std::vector<uint32_t> CommunityIndex; // The index in the vector is the node id and the value is the community id.


// !!! Assumes the iput ids start counting at 0 !!! 
EdgeList parseEdgeList(const std::string& path, uint32_t& nodeNum) {
	std::ifstream stream(path);

	if (!stream) {
		throw std::runtime_error("Failed to open file!");
		return EdgeList();
	}

	nodeNum = 0; // Assumption that indices start at 0.
	EdgeList edgeList{};

	std::string line;

	while (std::getline(stream, line)) {
		if (line.empty()) {
			continue;
		}

		std::stringstream stringStream(line);

		int a, b;

		if (stringStream >> a >> b) { // If there are two whitespace separated integer on this line.
			if (a > nodeNum) { // Take max.
				nodeNum = a;
			}

			if (b > nodeNum) { // Take max.
				nodeNum = b;
			}

			// IF zero based -> (a,b), IF one based -> (a-1,b-1)
			edgeList.push_back(std::make_pair(a, b));
		}
	}

	if (edgeList.empty()) {
		return edgeList;
	}

	++nodeNum;

	return edgeList;
}

AdjacencyList edgeListToAdjacencyList(const EdgeList& edgeList, const int nodeNum) {
	AdjacencyList adjacencyList(nodeNum);

	// Creaty empty adjacency lists.
	for (int i = 0; i < adjacencyList.size(); ++i) {
		adjacencyList[i] = std::vector<uint32_t>();
	}

	for (const auto& edge : edgeList) {
		adjacencyList[edge.first].push_back(edge.second);
		adjacencyList[edge.second].push_back(edge.first);
	}

	return adjacencyList;
}

template <typename T>
inline bool contains(const T* data, uint32_t size, const T& value) {
	for (int i = 0; i < size; ++i) {
		if (data[i] == value) {
			return true;
		}
	}

	return false;
}
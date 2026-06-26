#pragma once

#include <string>
#include <cstdlib>
#include <chrono>
#include <point.hpp>
#include <quadtree.hpp>
#include <projection_simulation.hpp>
#include <common.hpp>
#include <dbscan.hpp>

constexpr uint8_t BLACKHOLE_DIMENSIONS = 3;
constexpr uint32_t BLACKHOLE_MIN_PTS = 7;
constexpr float BLACKHOLE_PRUNING_FRACTION = 0.0f;

const std::string inputPath = "..\\input\\10000.dat";
const std::string outputPath = "..\\output\\result.dat";

int main(){
	using Clock = std::chrono::high_resolution_clock;

	auto start = Clock::now();

	srand((unsigned int) time( NULL));

	std::cout << "BlackHole Parameters: " << '\n';

	std::cout<<"Dimensions: "<< BLACKHOLE_DIMENSIONS << '\n';
	
	std::cout<<"MinPts: = "<< BLACKHOLE_MIN_PTS << '\n';

	std::cout<<"Pruning Factor: "<< BLACKHOLE_PRUNING_FRACTION << '\n';

	uint32_t nodeNum = 0;
	EdgeList edgeList = parseEdgeList(inputPath, nodeNum);

	AdjacencyList adjacencyList = edgeListToAdjacencyList(edgeList, nodeNum);

	ProjectionSimulation<BLACKHOLE_DIMENSIONS> projectionSimulation(adjacencyList, false);

	projectionSimulation.simulateAllIterations();

	// DBSCAN
	Point<BLACKHOLE_DIMENSIONS>* points = projectionSimulation.getResult();
	int* clusters = new	int[nodeNum]; // The resulting communities.

	auto startDbscan = Clock::now();

	dbscan(points, clusters, nodeNum, BLACKHOLE_MIN_PTS, estimateEpsilon(points, nodeNum, BLACKHOLE_MIN_PTS, BLACKHOLE_PRUNING_FRACTION));
		
	auto end = Clock::now();
	auto projectionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(startDbscan - start);
	auto dbscanDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end - startDbscan);
	auto bhDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	
	std::cout << "BlackHole done, it took: " << bhDuration.count() << "ms     Time spend on projecting: " << projectionDuration.count() << "ms     Time spent on clustering: " << dbscanDuration.count() << "ms\n";

	return 0;
}

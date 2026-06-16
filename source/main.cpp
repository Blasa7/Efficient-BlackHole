#pragma once

#include "edgeReader.hpp"
#include "nodeCollection.hpp"
#include "ClusterPlay.hpp"
#include "DBscanPlay.hpp"
#include <string>
#include <cstdlib>
#include <chrono>
#include <point.hpp>
#include <quadtree.hpp>
#include <projection_simulation.hpp>
#include <common.hpp>

constexpr uint32_t BLACKHOLE_DIMENSIONS = 3;
constexpr float BLACKHOLE_ALPHA = 0.0f;
constexpr uint32_t BLACKHOLE_MIN_PTS = 7;
constexpr float BLACKHOLE_PRUNING_FRACTION = 0.0f;

const std::string inputPath = "..\\input\\10000.dat";
const std::string outputPath = "..\\output\\result.dat";

int main(){
	Point<1> t = Point<1>();
	//Quadtree<1, 10> t1;

	using Clock = std::chrono::high_resolution_clock;

	auto start = Clock::now();

	srand((unsigned int) time( NULL));

	std::cout << "BlackHole Parameters: " << '\n';

	std::cout<<"Alpha: "<< BLACKHOLE_ALPHA<< '\n';

	std::cout<<"Dimensions: "<< BLACKHOLE_DIMENSIONS << '\n';
	
	std::cout<<"MinPts: = "<< BLACKHOLE_MIN_PTS << '\n';

	std::cout<<"Pruning Factor: "<< BLACKHOLE_PRUNING_FRACTION << '\n';

	//ClusterPlay<BLACKHOLE_DIMENSIONS> cp;

	//cp.play(inputPath.c_str(), BLACKHOLE_ALPHA, outputPath.c_str());

	uint32_t nodeNum = 0;
	EdgeList edgeList = parseEdgeList(inputPath, nodeNum);

	AdjacencyList adjacencyList = edgeListToAdjacencyList(edgeList, nodeNum);

	ProjectionSimulation<BLACKHOLE_DIMENSIONS> projectionSimulation(adjacencyList, false);

	projectionSimulation.simulateAllIterations();

	// DBSCAN
	DBscanPlay<BLACKHOLE_DIMENSIONS>::dbscanCalculator(projectionSimulation.getResult(), nodeNum, BLACKHOLE_MIN_PTS, BLACKHOLE_PRUNING_FRACTION);
	
	auto end = Clock::now();
	auto bhduration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	
	std::cout << "BlackHole done, it took: " << bhduration.count() << "ms! \n";

	return 0;
}

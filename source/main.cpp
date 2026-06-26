#pragma once

#include "edgeReader.hpp"
#include "nodeCollection.hpp"
#include "ClusterPlay.hpp"
#include "DBscanPlay.hpp"
#include <string>
#include <cstdlib>
#include <chrono>

constexpr uint32_t BLACKHOLE_DIMENSIONS = 3;
constexpr float BLACKHOLE_ALPHA = 0.05f;
constexpr uint32_t BLACKHOLE_MIN_PTS = 7;
constexpr float BLACKHOLE_PRUNING_FRACTION = 0.0f;

const std::string inputPath = "..\\input\\100000.dat";
const std::string outputPath = "..\\output\\result.dat";

int main(){
	using Clock = std::chrono::high_resolution_clock;

	auto start = Clock::now();

	srand((unsigned int) time( NULL));

	std::cout << "BlackHole Parameters: " << '\n';

	std::cout<<"Alpha: "<< BLACKHOLE_ALPHA<< '\n';

	std::cout<<"Dimensions: "<< BLACKHOLE_DIMENSIONS << '\n';
	
	std::cout<<"MinPts: = "<< BLACKHOLE_MIN_PTS << '\n';

	std::cout<<"Pruning Factor: "<< BLACKHOLE_PRUNING_FRACTION << '\n';

	ClusterPlay<BLACKHOLE_DIMENSIONS> cp;

	cp.play(inputPath.c_str(), BLACKHOLE_ALPHA, outputPath.c_str());

	auto startDbscan = Clock::now();

	// DBSCAN
	DBscanPlay<BLACKHOLE_DIMENSIONS>::dbscanCalculator(outputPath.c_str(), BLACKHOLE_MIN_PTS, BLACKHOLE_PRUNING_FRACTION);
	
	auto end = Clock::now();

	auto projectionDuration = std::chrono::duration_cast<std::chrono::milliseconds>(startDbscan - start);
	auto clusteringDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end - startDbscan);
	auto totalDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	
	std::cout << "BlackHole done, it took: " << totalDuration.count() << "ms     Time spend on projecting: " << projectionDuration.count() << "ms     Time spent on clustering: " << clusteringDuration.count() << "ms\n";

	return 0;
}

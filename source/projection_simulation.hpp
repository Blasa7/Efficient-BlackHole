#pragma once

#include <cstdint>
#include <vector>
#include <iostream>

#include <common.hpp>
#include <point.hpp>
#include <random.hpp>
#include <quadtree.hpp>

#define MAX_QUADTREE_DEPTH 18

#undef min
#undef max

// Class to perform modfidied Barnes-Hutt simulation to project a graph into Euclidian space for the BlackHole algorithm.
template <uint32_t dimensions>
class ProjectionSimulation {
public:
	ProjectionSimulation(const AdjacencyList& adjacencyList, bool storeIntermediateResults) : adjacencyList(adjacencyList), storeIntermediateResults(storeIntermediateResults) {		
		// First uniformly distribute the points.
		Point<dimensions>* points = static_cast<Point<dimensions>*>(std::malloc(adjacencyList.size() * sizeof(Point<dimensions>)));;// new Point<dimensions>[adjacencyList.size()];

		for (int i = 0; i < adjacencyList.size(); ++i) {
			for (int j = 0; j < dimensions; j++) {
				points[i][j] = random::nextf() - 0.5f; // Uniform random value in range [-0.5, 0.5]
			}
		}

		results.push_back(points);

		// Set the repulsion factor. Could pass the total degree if it is known beforehand.
		uint32_t totalDegree = 0;
		for (int i = 0; i < adjacencyList.size(); ++i) {
			totalDegree += adjacencyList[i].size();
		}

		repulsionFactor = 1.0f / totalDegree * std::powf(totalDegree, 0.5f * (attractionExponent - repulsionExponent));

		// Must find the minCoordinate and maxCoordinate for every iteration as we rebuild the quadtree. Could also search for it in the octree from the previuous iteration.
		minCoordinate = points[0];
		maxCoordinate = points[0];

		for (int i = 0; i < adjacencyList.size(); ++i) {
			for (int j = 0; j < dimensions; ++j) {
				if (points[i][j] < minCoordinate[j]) {
					minCoordinate[j] = points[i][j];
				}

				if (points[i][j] > maxCoordinate[j]) {
					maxCoordinate[j] = points[i][j];
				}
			}
		}

		// Provide padding for the points to move around.
		for (int i = 0; i < dimensions; ++i) {
			float padding = (maxCoordinate[i] - minCoordinate[i]) / 2;

			maxCoordinate[i] += padding;
			minCoordinate[i] -= padding;
		}

		std::cout << "initial maximum coordinate x: " << std::to_string(maxCoordinate[0]) << '\n';

		// Build the quadtree.
		quadtree = new Quadtree<dimensions, MAX_QUADTREE_DEPTH>(minCoordinate, maxCoordinate, 0, nullptr);

		// Calculate initial energy.
		handles.reserve(adjacencyList.size());

		for (int i = 0; i < adjacencyList.size(); ++i) {
			handles.push_back(quadtree->addPoint(points[i], adjacencyList[i].size()));
		}

		for (int i = 0; i < adjacencyList.size(); ++i) {
			totalEnergy += quadtree->calculateRepulsiveEnergy(points[i], repulsionFactor, adjacencyList[i].size(), repulsionExponent) + calculateAttractiveEnergy(i);
		}

		std::cout << "Initial Energy: " << totalEnergy << '\n';
	}

	~ProjectionSimulation() {
		for (int i = 0; i < results.size(); ++i) {
			std::free(results[i]);
		}
	}

	void simulateAllIterations() {
		for (int i = iteration; i < maxIterations; ++i) {
			simulateIteration();

			std::cout << "Iteration: " << i <<"     Energy: " << totalEnergy << "     Attraction: " << attractionExponent << "     Repulsion : " << repulsionExponent << '\n';
		}
	}

	void simulateIteration() {
		totalEnergy = 0.0f;

		Point<dimensions>* points = nullptr;

		if (storeIntermediateResults) {
			points = static_cast<Point<dimensions>*>(std::malloc(adjacencyList.size() * sizeof(Point<dimensions>)));
			
			std::memcpy(points, results.back(), adjacencyList.size() * sizeof(Point<dimensions>));

			results.push_back(points);
		}
		else {
			points = results.back();
		}

		// Must find the minCoordinate and maxCoordinate for every iteration as we rebuild the quadtree. Could also search for it in the octree from the previuous iteration.
		Point<dimensions> newMinCoordinate = points[0];
		Point<dimensions> newMaxCoordinate = points[0];
		
		bool remakeQuadTree = false;

		for (int i = 0; i < adjacencyList.size(); ++i) {
			for (int j = 0; j < dimensions; ++j) {
				if (points[i][j] < minCoordinate[j]) {
					minCoordinate[j] = points[i][j];
					remakeQuadTree = true;
				}

				if (points[i][j] > maxCoordinate[j]) {
					maxCoordinate[j] = points[i][j];
					remakeQuadTree = true;
				}
			}
		}


		// Rebuild the quadtree. Performance cost is almost nonexistent as it generally 
		// gives more balanced trees and the construction cost is not too large.

		// Provide padding for the points to move around.
		for (int i = 0; i < dimensions; ++i) {
			float padding = (maxCoordinate[i] - minCoordinate[i]) / 2;

			newMaxCoordinate[i] += padding;
			newMinCoordinate[i] -= padding;
		}

		delete quadtree;

		quadtree = new Quadtree<dimensions, MAX_QUADTREE_DEPTH>(newMinCoordinate, newMaxCoordinate, 0, nullptr);

		minCoordinate = newMinCoordinate;
		maxCoordinate = newMaxCoordinate;

		for (int i = 0; i < adjacencyList.size(); ++i) {
			handles[i] = quadtree->addPoint(points[i], adjacencyList[i].size());
		}
		
		adjustExponents();

		// Simulate the movement of every point in series.
		for (int i = 0; i < adjacencyList.size(); ++i) {
			Point<dimensions>& point = points[i];

			// First calculate the direction the forces want to move the point in.
			Point<dimensions> direction;

			// Repulsive direction contribution.
			direction.add(quadtree->calculateRepulsiveDirection(point, repulsionFactor, adjacencyList[i].size(), repulsionExponent));

			// Attractive direction contribution.
			for (int j = 0; j < adjacencyList[i].size(); ++j) {
				Point<dimensions>& neighbor = points[adjacencyList[i][j]];

				float squaredDistance = Point<dimensions>::squaredDistance(point, neighbor);

				if (squaredDistance == 0.0f) {
					continue;
				}

				Point<dimensions> neighborDirection;

				for (int k = 0; k < dimensions; ++k) {
					neighborDirection[k] = (neighbor[k] - point[k]);
				}

				// We use the squared distance instead of the normal distance and thus add a factor 0.5f to the exponent.
				neighborDirection.scale(std::powf(squaredDistance, 0.5f * attractionExponent - 1.0f));

				direction.add(neighborDirection);
			}

			direction.scale(1.0f / 32.0f);

			float weighting = adjacencyList[i].size();

			int bestScale = 0;
			float bestEnergy = quadtree->calculateRepulsiveEnergy(point, repulsionFactor, weighting, repulsionExponent) + calculateAttractiveEnergy(i);

			Point<dimensions> oldPosition = point;

			// Decrease scale to find best candidate.
			// 1. remove point
			// 2. compute new position
			// 3. add point back
			// 4. calculate and compare new energy
			for (int scale = 32; scale >= 1 && (bestScale == 0 || bestScale / 2 == scale); scale /= 2) {
				handles[i].removePoint(point, weighting);

				for (int j = 0; j < dimensions; ++j) {
					point[j] = oldPosition[j] + direction[j] * scale;
				}

				handles[i] = quadtree->addPoint(point, weighting);

				// Fairly expensive energy calculation call, could look into incremental energy delta computation.
				float newEnergy = quadtree->calculateRepulsiveEnergy(point, repulsionFactor, weighting, repulsionExponent) + calculateAttractiveEnergy(i);

				if (newEnergy < bestEnergy) {
					bestEnergy = newEnergy;
					bestScale = scale;
				}
			}

			// Increase scale to find best candidate if the best scale was 32.
			// 1. remove point
			// 2. compute new position
			// 3. add point back
			// 4. calculate and compare new energy
			for (int scale = 64; scale <= 128 && bestScale == scale / 2; scale *= 2) {
				handles[i].removePoint(point, weighting);

				for (int j = 0; j < dimensions; ++j) {
					point[j] = oldPosition[j] + direction[j] * scale;
				}

				handles[i] = quadtree->addPoint(point, weighting);

				// Fairly expensive energy calculation call, could look into incremental energy delta computation.
				float newEnergy = quadtree->calculateRepulsiveEnergy(point, repulsionFactor, weighting, repulsionExponent) + calculateAttractiveEnergy(i);

				if (newEnergy < bestEnergy) {
					bestEnergy = newEnergy;
					bestScale = scale;
				}
			}

			// Remove and add back with the best scale.
			handles[i].removePoint(point, weighting);

			for (int j = 0; j < dimensions; ++j) {
				point[j] = oldPosition[j] + direction[j] * bestScale;
			}

			handles[i] = quadtree->addPoint(point, weighting);

			totalEnergy += bestEnergy;
		}

		++iteration;
	}

	// Returns all stored projections.
	std::vector<Point<dimensions>*> getResults() {
		return results;
	}

	// Returns the last stored projection.
	Point<dimensions>* getResult() {
		return results.back();
	}
private:
	const AdjacencyList& adjacencyList;

	bool storeIntermediateResults;
	
	std::vector<Point<dimensions>*> results;

	float repulsionFactor; // Set in constructor.

	// BlackHole constants.
	float attractionExponent = 0.00f;
	float repulsionExponent = 0.0f;

	float attractionExponentFinal = 0.00f;
	float repulsionExponentFinal = 0.0f;

	uint32_t maxIterations = 100;
	uint32_t iteration = 0;

	float totalEnergy = 0.0f;

	Point<dimensions> minCoordinate;
	Point<dimensions> maxCoordinate;

	// Reuse between iterations.
	Quadtree<dimensions, MAX_QUADTREE_DEPTH>* quadtree = nullptr;
	std::vector<QuadtreePointHandle<dimensions, MAX_QUADTREE_DEPTH>> handles;


	

	void adjustExponents() {
		uint32_t i = iteration + 1;

		if (maxIterations >= 50 && repulsionExponentFinal < 1.0f) {
			attractionExponent = attractionExponentFinal;
			repulsionExponent = repulsionExponentFinal;

			if (i <= 0.6f * maxIterations) {
				attractionExponent = attractionExponent + 1.1f * (1.0f - repulsionExponentFinal);
				repulsionExponent = repulsionExponent + 0.9f * (1.0f - repulsionExponentFinal);
			}
			else if (i <= 0.9f * maxIterations) {
				attractionExponent = attractionExponent + 1.1f * (1.0f - repulsionExponentFinal) * (0.9f - (i / (float)maxIterations)) / 0.3f;
				repulsionExponent = repulsionExponent + 0.9f * (1.0f - repulsionExponentFinal) * (0.9f - (i / (float)maxIterations)) / 0.3f;
			}
		}
	}

	float calculateAttractiveEnergy(const uint32_t index) {
		const Point<dimensions>& point = results.back()[index];

		float energy = 0.0f;

		// Attractive direction contribution.
		for (int j = 0; j < adjacencyList[index].size(); ++j) {
			const Point<dimensions>& neighbor = results.back()[adjacencyList[index][j]];

			float squaredDistance = Point<dimensions>::squaredDistance(point, neighbor);

			if (attractionExponent == 0.0f) {
				// We use the squared distance instead of the normal distance and thus add a factor 0.5f.
				energy += 0.5f * std::logf(squaredDistance);
			}
			else {
				// We use the squared distance instead of the normal distance and thus add a factor 0.5f.
				energy += std::powf(squaredDistance, 0.5f * attractionExponent) / attractionExponent;
			}
		}

		return energy;
	}
};
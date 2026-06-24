#pragma once

#include <point.hpp>
#include <kd_tree.hpp>

#include <algorithm>
#include <limits>
#include <vector>

// Will modify points.
template <uint32_t dimensions>
void dbscan(const Point<dimensions>* inputPoints, int* clusters, int pointNum, int minPts, float epsilon) {
	float epsilonSquared = epsilon * epsilon;

	// Necessary because the kdtree will reorder the points and we need to retrieve the original index.
	PointWithID<dimensions>* points = static_cast<PointWithID<dimensions>*>(malloc(pointNum * sizeof(PointWithID<dimensions>)));

	for (int i = 0; i < pointNum; ++i) {
		clusters[i] = -1; // Mark unvisited.

		new (&points[i]) PointWithID<dimensions>(inputPoints[i], i);
	}

	int nextClusterID = 0;

	KDTree<dimensions> pointLookup(points, 0, pointNum - 1, 0);

	for (int i = 0; i < pointNum; ++i) {
		PointWithID<dimensions> point = points[i];

		// If unvisited.
		if (clusters[point.id] == -1) {
			std::vector<PointWithID<dimensions>> neighbors;
			pointLookup.rangeSearchSquaredDistance(point, epsilonSquared, neighbors);

			// Cluster size too small, label as noise.
			if (neighbors.size() < minPts) {
				clusters[point.id] = -2;
				continue;
			}

			// A new cluster is found.
			clusters[point.id] = nextClusterID;

			while (neighbors.size() > 0) {
				PointWithID<dimensions> neighbor = neighbors.back();
				neighbors.pop_back();
				
				// If its labelled as noise then assign it to the cluster and dont search its neighbors.
				if (clusters[neighbor.id] == -2) {
					clusters[neighbor.id] = nextClusterID;
					continue;
				}

				// If it already assigned to a different cluster (border point) then do nothing.
				if (clusters[neighbor.id] >= 0) {
					continue;
				}

				// Otherwise add it to the cluster and its neighbors to the unexplored neighbor queue.
				clusters[neighbor.id] = nextClusterID;

				std::vector<PointWithID<dimensions>> nextNeighbors;
				pointLookup.rangeSearchSquaredDistance(neighbor, epsilonSquared, nextNeighbors);

				neighbors.insert(neighbors.end(), nextNeighbors.begin(), nextNeighbors.end());
			}

			++nextClusterID;
		}
	}

	// Any remaining noise points are assigned into one big community just like the original paper. 
	for (int i = 0; i < pointNum; ++i) {
		if (clusters[i] < 0) {
			clusters[i] = nextClusterID;
		}
	}

	std::free(points);
}
template<uint32_t dimensions>
float estimateEpsilon(Point<dimensions>* points, int pointNum, int minPts, float pruningFraction) {
	float* x = new float[pointNum];
	float* y = new float[pointNum];

	// Use x for distances temporarily.
	for (int i = 0; i < pointNum; ++i) {
		for (int j = 0; j < pointNum; ++j) {
			x[j] = Point<dimensions>::distance(points[i], points[j]);// calcDistance<dimensions>(points[i], points[j]);
		}

		std::nth_element(x, x + minPts - 1, x + pointNum);
		y[i] = x[minPts - 1];
	}

	// Sort ascending instead.
	std::sort(y, y + pointNum);

	// Then trunctuate at the end. We dont actually need to zero the values.
	int remainingPointNum = (int)(pointNum * (1.0f - pruningFraction));

	float* maxDistances = new float[remainingPointNum];

	// Save only the used values.
	for (int i = 0; i < remainingPointNum; ++i) {
		maxDistances[i] = y[i];
	}

	float minY = std::numeric_limits<float>::max();
	float maxY = std::numeric_limits<float>::min();

	// Flipped because the order is flipped.
	for (int i = 0; i < remainingPointNum; ++i) {
		x[i] = remainingPointNum - i - 1;

		if (y[i] < minY) {
			minY = y[i];
		}

		if (y[i] > maxY) {
			maxY = y[i];
		}
	}

	for (int i = 0; i < remainingPointNum; ++i) {
		x[i] = x[i] / pointNum;
		y[i] = (y[i] - minY) / (maxY - minY);
	}

	for (int i = 0; i < remainingPointNum; ++i) {
		x[i] = cosf(-PI / 4.0f) * x[i] + sinf(-PI / 4.0f) * (y[i] - 1.0f);
		y[i] = -sin(-PI / 4.0f) * x[i] + cos(-PI / 4.0f) * (y[i] - 1.0f);
	}

	minY = std::numeric_limits<float>::max();
	int epsIndex = -1;

	for (int i = 0; i < remainingPointNum; ++i) {
		if (y[i] < minY) {
			minY = y[i];
			epsIndex = i;
		}
	}

	float eps = maxDistances[epsIndex];

	delete[] maxDistances;
	delete[] x;
	delete[] y;

	return eps;
}
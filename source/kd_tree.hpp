#pragma once

#include <dbscan.hpp>
#include <random.hpp>
#include <point.hpp>

#include <algorithm>
#include <vector>

template <uint32_t dimensions>
class PointWithID : public Point<dimensions> {
public:
	uint32_t id;

	PointWithID(Point<dimensions> point, uint32_t id) : Point<dimensions>(point.data()), id(id) {		
		
	}

	PointWithID() : Point<dimensions>(), id(0) {

	}
};

template<uint32_t dimensions>
class KDTree {
public:
	KDTree(PointWithID<dimensions>* points, const int startIndex, const int endIndex, int depth) {
		// No points
		if (startIndex == endIndex) { 
			left = nullptr;
			right = nullptr;

			point = points[startIndex];

			return;
		} 

		// This is not a leaf node, split a certain dimension based on depth.
		const int medianIndex = startIndex + ((endIndex - startIndex) >> 1);
		const int splitDimension = depth % dimensions;

		point = quickSelect(points, startIndex, endIndex, medianIndex, splitDimension);
		
		++depth;

		// Create left child if possible.
		if (medianIndex == startIndex) {
			left = nullptr;
		}
		else {
			left = new KDTree(points, startIndex, medianIndex - 1, depth);
		}

		if (medianIndex == endIndex) {
			right = nullptr;
		}
		else
		{
			right = new KDTree(points, medianIndex + 1, endIndex, depth);
		}
	}

	~KDTree() {
		delete left;
		delete right;
	}

	void rangeSearchSquaredDistance(PointWithID<dimensions> centre, float maximumSquaredDistance, std::vector<PointWithID<dimensions>>& result, bool* inQueue) {
		return rangeSearchSquaredDistance(centre, maximumSquaredDistance, result, 0, inQueue);
	}

private:
	KDTree* left;
	KDTree* right;
	
	PointWithID<dimensions> point; // Median point or the leaf point. Depending on if this is a node or a leaf.

	// Quick select for expected linear time median selection.
	inline PointWithID<dimensions> quickSelect(PointWithID<dimensions> points[], const int startIndex, const int endIndex, const int medianIndex, const int splitDimension) {
		if (startIndex == endIndex) {
			return points[startIndex];
		}

		int pivotIndex = startIndex + random::range(0, endIndex - startIndex);
		
		pivotIndex = partition(points, startIndex, endIndex, splitDimension, pivotIndex);

		if (pivotIndex == medianIndex) {
			return points[medianIndex];
		}
		else if (medianIndex < pivotIndex) {
			return quickSelect(points, startIndex, pivotIndex - 1, medianIndex, splitDimension);
		}
		else {
			return quickSelect(points, pivotIndex + 1, endIndex, medianIndex, splitDimension);
		}
	}

	// Partitions a segment of point for the given dimension and pivot point, then returns the final pivot index.
	inline int partition(PointWithID<dimensions> points[], const int startIndex, const int endIndex, const int splitDimension,  const int pivotIndex) {
		const float pivotValue = points[pivotIndex][splitDimension];

		std::swap(points[endIndex], points[pivotIndex]);

		int store = startIndex;

		for (int i = startIndex; i < endIndex; ++i) {
			if (points[i][splitDimension] < pivotValue) {
				std::swap(points[store], points[i]);
				++store;
			}
		}

		std::swap(points[store], points[endIndex]);

		return store;
	}

	inline void rangeSearchSquaredDistance(const PointWithID<dimensions>& centre, float maximumSquaredDistance, std::vector<PointWithID<dimensions>>& result, int depth, bool* inQueue) {
		if (PointWithID<dimensions>::squaredDistance(centre, point) <= maximumSquaredDistance && !inQueue[point.id]) {
			result.push_back(point);
			inQueue[point.id] = true;
		}

		int splitDimension = depth % dimensions;

		float diff = centre[splitDimension] - point[splitDimension];
		
		++depth;

		// Negative so centre is on the left of the plane and the left branch is the near branch.
		if (diff < 0) {
			if (left != nullptr) {
				left->rangeSearchSquaredDistance(centre, maximumSquaredDistance, result, depth, inQueue);
			}

			if (diff * diff <= maximumSquaredDistance && right != nullptr) { // Must also search right of the plane.
				right->rangeSearchSquaredDistance(centre, maximumSquaredDistance, result, depth, inQueue);
			}
		} 
		else { // Positive so centre is on the right side of the plane and the right branch is the near branch.
			if (right != nullptr) {
				right->rangeSearchSquaredDistance(centre, maximumSquaredDistance, result, depth, inQueue);

			}

			if (diff * diff <= maximumSquaredDistance && left != nullptr) { // Must also search right of the plane.
				left->rangeSearchSquaredDistance(centre, maximumSquaredDistance, result, depth, inQueue);
			}
		}
	}
};
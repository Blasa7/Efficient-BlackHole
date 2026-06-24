#pragma once

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <common.hpp>
#include <point.hpp>

template <uint8_t dimensions, uint8_t maxDepth>
class QuadtreePointHandle;

// Quadtree for arbitrary dimensions.
template <uint8_t dimensions, uint8_t maxDepth>
class Quadtree {
public:
	Quadtree(const Point<dimensions>& minCoordinates, const Point<dimensions>& maxCoordinates, uint8_t depth, Quadtree<dimensions, maxDepth>* parent) : minCoordinates(minCoordinates), maxCoordinates(maxCoordinates), depth(depth), parent(parent) {
		// Only nodes that are not at max depth need space allocated for child nodes.
		if (depth != maxDepth) {
			// We are only allocating the space and not actually constructor the children for two reasons. 
			// 1. Many children may never actually get used. 
			// 2. If children get constructed immediately it would cascade down to a full tree containing childCount ^ maxDepth which is way more nodes than likely needed. 
			children = static_cast<Quadtree<dimensions, maxDepth>*>(std::malloc(childCount * sizeof(Quadtree<dimensions, maxDepth>)));

			// Initially no children exist and the mapping are cleared.
			for (int i = 0; i < childCount; ++i) {
				childMapping[i] = -1;
			}
		}

		float maxWidth = calculateMaxWidth();
		maxWidthSquared = maxWidth * maxWidth;
	}

	~Quadtree() {
		if (children != nullptr) {
			for (int i = 0; i < nextChildIndex; ++i) {
				children[i].~Quadtree();
			}

			std::free(children);
		}
	}

	QuadtreePointHandle<dimensions, maxDepth> addPoint(const Point<dimensions>& point, const float weighting) {
		QuadtreePointHandle<dimensions, maxDepth> result;

		// If nothing has been stored here or the maximum depth has been reached then this point must be stored here.
		if (depth == maxDepth || (weight <= 0.0f && isLeaf)) {
			result = QuadtreePointHandle<dimensions, maxDepth>(this);

			++numberOfBranchLeaves;
		}
		// Not maximum depth and no space to store the point. Thus we must split this node and store both its currently stored node and the new node in the children.
		else {

			// Move the stored node down into some child if it has not been done yet. Note that this invalidates the old handle somewhat and requires searching children again to find it.
			if (isLeaf) {
				int childIndexOld = 0;
				Point<dimensions> childMinCoordinatesOld;
				Point<dimensions> childMaxCoordinatesOld;

				for (int i = 0; i < dimensions; ++i) {
					float centre = (minCoordinates[i] + maxCoordinates[i]) / 2.0f;

					if (centreOfMass[i] > centre) { // TODO: store cental point instead.
						childIndexOld += 1 << i;
						childMinCoordinatesOld[i] = centre;
						childMaxCoordinatesOld[i] = maxCoordinates[i];
					}
					else {
						childMinCoordinatesOld[i] = minCoordinates[i];
						childMaxCoordinatesOld[i] = centre;
					}
				}

				// Create the child Quadtrees if they dont exist yet.
				if (childMapping[childIndexOld] == -1) {
					new (&children[nextChildIndex]) Quadtree<dimensions, maxDepth>(childMinCoordinatesOld, childMaxCoordinatesOld, depth + 1, this);
					childMapping[childIndexOld] = nextChildIndex;
					++nextChildIndex;
				}

				// Insert the points into the children and store the handle for efficient deletion.
				pointHandle = children[childMapping[childIndexOld]].addPoint(centreOfMass, weight);

				isLeaf = false;
			}

			// Now do the same but for the new point.

			int childIndexNew = 0;
			Point<dimensions> childMinCoordinatesNew;
			Point<dimensions> childMaxCoordinatesNew;

			for (int i = 0; i < dimensions; ++i) {
				float centre = (minCoordinates[i] + maxCoordinates[i]) / 2.0f;

				if (point[i] > centre) { // TODO: store cental point instead.
					childIndexNew += 1 << i;
					childMinCoordinatesNew[i] = centre;
					childMaxCoordinatesNew[i] = maxCoordinates[i];
				}
				else {
					childMinCoordinatesNew[i] = minCoordinates[i];
					childMaxCoordinatesNew[i] = centre;
				}
			}

			// Create the child Quadtrees if they dont exist yet.
			if (childMapping[childIndexNew] == -1) {
				new (&children[nextChildIndex]) Quadtree<dimensions, maxDepth>(childMinCoordinatesNew, childMaxCoordinatesNew, depth + 1, this);
				childMapping[childIndexNew] = nextChildIndex;
				++nextChildIndex;
			}

			// Insert the points into the children.
			result = children[childMapping[childIndexNew]].addPoint(point, weighting);
		}

		// Update the centre of mass. Note that if this first point is added the centre of mass is equal to that points position.
		// This is also used when splitting a leaf node that is not at maximum depth.
		for (int i = 0; i < dimensions; ++i) {
			centreOfMass[i] = (weight * centreOfMass[i] + weighting * point[i]) / (weight + weighting); // TODO, check if it can be simplified.
		}

		weight += weighting;

		return result;
	}

	float calculateRepulsiveEnergy(const Point<dimensions>& point, const float repulsionFactor, const float pointWeight, const float repulsionExponent) {
		return -repulsionFactor * pointWeight * calculateUnscaledRepulsiveEnergy(point, repulsionExponent);
	}

	Point<dimensions> calculateRepulsiveDirection(const Point<dimensions>& point, const float repulsionFactor, const float pointWeight, const float repulsionExponent) {
		Point<dimensions> direction = Point<dimensions>();

		calculateUnscaledRepulsiveDirection(point, repulsionExponent, direction);

		direction.scale(repulsionFactor * pointWeight);

		return direction;
	}


private:
	static constexpr uint8_t childCount = 1 << dimensions;

	Quadtree<dimensions, maxDepth>* parent;

	Quadtree<dimensions, maxDepth>* children = nullptr;
	int8_t childMapping[childCount];

	// Used to reduce recursion when the number of leaves in this branch equals one. As it means all points in this branch are in a single leaf.
	uint32_t numberOfBranchLeaves = 0;

	QuadtreePointHandle<dimensions, maxDepth> pointHandle; // If this node is split then the point is move to some child. This pointer points to that child.

	uint8_t nextChildIndex = 0;
	uint8_t depth;
	bool isLeaf = true;

	Point<dimensions> minCoordinates;
	Point<dimensions> maxCoordinates;
	float maxWidthSquared;

	Point<dimensions> centreOfMass = Point<dimensions>(); // Zero constructor.
	float weight = 0.0;

	// Used as entry by the handle to remove a point.
	void removePoint(const Point<dimensions>& point, const float weighting) {
		// It possible tht the point for the handle is no longer stored here and has been moved downwards to some child on splitting. 
		// Thus we check for this and search the children if necesary. Note that the tree is guaranteed to be relatively shallow by the maxDepth parameter
		// meaning it should remain fast even in the worst cases. Additionally the majority of points will be stored at the bottom of the tree and never split.
		if (!pointHandle.isEmpty()) {

			pointHandle.removePoint(point, weighting);

			return; // Removal will be handled by upwards propagation.
		}

		uint32_t leavesRemoved;

		if (weight - weighting <= 0.0f) {
			leavesRemoved = 1;
		}

		removePointUpwards(point, weighting, leavesRemoved); // Propagate removal upwards.
	}

	// Used to propagate upwards point removal.
	void removePointUpwards(const Point<dimensions>& point, const float weighting, const uint32_t leavesRemoved) {
		// If we are trying to remove the last point, reset the centre of mass to 0 to avoid division by 0 corrupting future reuse of the centre of mass.
		if (weight <= weighting) {
			centreOfMass = Point<dimensions>();
		}
		// Otherwise recalculate for removing the given point.
		else {
			for (int i = 0; i < dimensions; ++i) {
				centreOfMass[i] = (weight * centreOfMass[i] - weighting * point[i]) / (weight - weighting);
			}
		}

		weight -= weighting;
		numberOfBranchLeaves -= leavesRemoved;

		if (parent != nullptr) {
			parent->removePointUpwards(point, weighting, leavesRemoved);
		}
	}

	// To propagate leaf creation upwards.
	void leafCreated() {
		++numberOfBranchLeaves;

		if (parent != nullptr) {
			parent->leafCreated();
		}
	}

	// Calculates the repulsive energy omitting the constant multipliers.
	float calculateUnscaledRepulsiveEnergy(const Point<dimensions>& point, const float repulsiveExponent) const {
		float squaredDistance = Point<dimensions>::squaredDistance(point, centreOfMass);

		// A point may not repel itself.
		if (squaredDistance == 0.0f) {
			return 0.0f;
		}

		// If we can consider this branch as a single mass. Leaf nodes may also be considered as a single mass. 
		// Leaf nodes may either be at maximum depth or contain exactly one point.
		if (squaredDistance >= maxWidthSquared || isLeaf || numberOfBranchLeaves == 1) {
			if (repulsiveExponent == 0.0f) {
				// We use the squared distance instead of the normal distance and thus add a factor 0.5f.
				return weight * 0.5f * std::logf(squaredDistance);
			}
			else {
				// We use the squared distance instead of the normal distance and thus add a factor 0.5f.
				return weight * std::powf(squaredDistance, 0.5f * repulsiveExponent) / repulsiveExponent;
			}
		}

		// Otherwise this must not be a leaf node or be treade as a single mass and we must recurse into all child trees..
		float energy = 0.0f;

		for (int i = 0; i < nextChildIndex; ++i) {
			if (children[i].weight > 0.0f) {
				energy += children[i].calculateUnscaledRepulsiveEnergy(point, repulsiveExponent);
			}
		}

		return energy;
	}

	// Calculates the repulsive direction omitting the constant multipliers.
	void calculateUnscaledRepulsiveDirection(const Point<dimensions>& point, const float repulsiveExponent, Point<dimensions>& directionAccumulator) const {
		float squaredDistance = Point<dimensions>::squaredDistance(point, centreOfMass);

		if (squaredDistance == 0.0f) {
			return;
		}

		// If we can consider this branch as a single mass. Leaf nodes may also be considered as a single mass. 
		// Leaf nodes may either be at maximum depth or contain exactly one point.
		if (squaredDistance >= maxWidthSquared || isLeaf || numberOfBranchLeaves == 1) {
			// We use the squared distance instead of the normal distance and thus add a factor 0.5f to the exponent.
			float scale = weight * std::powf(squaredDistance, 0.5f * repulsiveExponent - 1.0f);

			for (int i = 0; i < dimensions; ++i) {
				directionAccumulator[i] += (point[i] - centreOfMass[i]) * scale;
			}

			return;
		}

		// Otherwise this must not be a leaf node or be treated as a single mass and we must recurse into all child trees..
		for (int i = 0; i < nextChildIndex; ++i) {
			if (children[i].weight > 0.0f) {
				children[i].calculateUnscaledRepulsiveDirection(point, repulsiveExponent, directionAccumulator);
			}
		}
	}

	float calculateMaxWidth() {
		float maxWidth = 0.0f;

		for (int i = 0; i < dimensions; ++i) {
			float width = maxCoordinates[i] - minCoordinates[i];
			if (width > maxWidth) {
				maxWidth = width;
			}
		}

		return maxWidth;
	}

	// Dont allow reassignment as we are manually managing some memory that would not get freed in such cases corrupting the heap.
	Quadtree& operator=(const Quadtree&) = delete;

	// Allow the handle class to remove points.
	friend class QuadtreePointHandle<dimensions, maxDepth>;
};


// Handle generated by Quadtree addPoint function. Allows for efficient point deletion. A point handle should never be used twice to remove a point.
template <uint8_t dimensions, uint8_t maxDepth>
class QuadtreePointHandle {
public:
	void removePoint(const Point<dimensions>& point, const float weighting) {
		container->removePoint(point, weighting);
	}

private:
	Quadtree<dimensions, maxDepth>* container;

	QuadtreePointHandle(Quadtree<dimensions, maxDepth>* container) : container(container) {

	}

	// Empty handle constructor.
	QuadtreePointHandle() : container(nullptr) {

	}

	bool isEmpty() {
		return container == nullptr;
	}

	friend Quadtree<dimensions, maxDepth>;
};
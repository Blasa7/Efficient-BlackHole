#pragma once

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <common.hpp>
#include <point.hpp>

template <uint32_t dimensions, uint32_t maxDepth>
class QuadtreePointHandle;

// Quadtree for arbitrary dimensions.
template <uint32_t dimensions, uint32_t maxDepth>
class Quadtree {
public:
	Quadtree(const Point<dimensions>& minCoordinates, const Point<dimensions>& maxCoordinates, uint32_t depth, Quadtree<dimensions, maxDepth>* parent) : minCoordinates(minCoordinates), maxCoordinates(maxCoordinates), depth(depth), parent(parent) {
		// At maximum depth we have points but we do not have any child Quadtrees.
		if (depth == maxDepth) {
			maxPointCount = 4; // Preallocate a little bit more space.
		}
		// We have child Quadtrees but no more than 1 point.
		else {
			// We are only allocating the space and not actually constructor the children for two reasons. 
			// 1. Many children may never actually get used. 
			// 2. If children get constructed immediately it would cascade down to a full tree containing childCount ^ maxDepth which is way more nodes than likely needed. 
			children = static_cast<Quadtree<dimensions, maxDepth>*>(std::malloc(childCount * sizeof(Quadtree<dimensions, maxDepth>)));
		}

		points = static_cast<Point<dimensions>*>(std::malloc(maxPointCount * sizeof(Point<dimensions>)));

		float maxWidth = calculateMaxWidth();
		maxWidthSquared = maxWidth * maxWidth;
	}

	~Quadtree() {
		std::free(points);

		if (children != nullptr) {
			for (int i = 0; i < childCount; ++i) {
				//delete children[i];
				if (childCreated[i]) {
					children[i].~Quadtree();
				}
			}

			std::free(children);
		}
	}

	QuadtreePointHandle<dimensions, maxDepth> addPoint(const Point<dimensions>& point, const float weighting) {
		QuadtreePointHandle<dimensions, maxDepth> result;

		// If nothing has been stored here or the maximum depth has been reached then this point must be stored here.
		if (depth == maxDepth || pointCount == 0 || freePointIndices.size() >= 1) {
			uint32_t id = addPointHere(point);

			result = QuadtreePointHandle<dimensions, maxDepth>(this, id);
		}
		// Not maximum depth and no space to store the point. Thus we must split this node and store both its currently stored node and the new node in the children.
		else {

			// Move the stored node down into some child if it has not been done yet. Note that this invalidates the old handle somewhat and requires searching children again to find it.
			if (maxPointCount != 0) {
				int childIndexOld = 0;
				Point<dimensions> childMinCoordinatesOld;
				Point<dimensions> childMaxCoordinatesOld;

				for (int i = 0; i < dimensions; ++i) {
					float centre = (minCoordinates[i] + maxCoordinates[i]) / 2.0f;

					if (points[0][i] > centre) { // TODO: store cental point instead.
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
				if (!childCreated[childIndexOld]) {
					new (&children[childIndexOld]) Quadtree<dimensions, maxDepth>(childMinCoordinatesOld, childMaxCoordinatesOld, depth + 1, this);
					childCreated[childIndexOld] = true;
				}

				// Insert the points into the children and store the handle for efficient deletion.
				pointHandle = children[childIndexOld].addPoint(points[0], weight);

				maxPointCount = 0; // Effectively marking this node as already split.
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
			if (!childCreated[childIndexNew]) {
				new (&children[childIndexNew]) Quadtree<dimensions, maxDepth>(childMinCoordinatesNew, childMaxCoordinatesNew, depth + 1, this);
				childCreated[childIndexNew] = true;
			}

			// Insert the points into the children.
			result = children[childIndexNew].addPoint(point, weighting);
		}

		// Update the centre of mass.
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
		Point<dimensions> direction = calculateUnscaledRepulsiveDirection(point, repulsionExponent);
		
		direction.scale(repulsionFactor * pointWeight);

		return direction;
	}


private:
	static constexpr uint32_t childCount = 1 << dimensions;

	Quadtree<dimensions, maxDepth>* parent;

	Quadtree<dimensions, maxDepth>* children = nullptr;
	bool childCreated[childCount] = {};

	Point<dimensions>* points = nullptr; // Only leaves will get values. Used with std::malloc, std::free and std::realloc, do not use delete[] or new[].

	QuadtreePointHandle<dimensions, maxDepth> pointHandle; // If this node is split then the point is move to some child. This pointer points to that child.

	uint32_t pointCount = 0;
	uint32_t maxPointCount = 1;

	uint32_t depth;

	Point<dimensions> minCoordinates;
	Point<dimensions> maxCoordinates;
	float maxWidthSquared;

	Point<dimensions> centreOfMass = Point<dimensions>(); // Zero constructor.
	float weight = 0.0;

	// Storing the indices of deleted points in the points array. Allows efficient reuse of memory.
	std::vector<uint32_t> freePointIndices;

	inline uint32_t addPointHere(const Point<dimensions>& point) {
		// If some space has been freed previusly we can reuse that space.
		if (freePointIndices.size() > 0) {
			uint32_t index = freePointIndices.back();
			points[index] = point;

			freePointIndices.pop_back();

			return index;
		}

		// Array is full, allocate more space.
		if (pointCount == maxPointCount) {
			maxPointCount <<= 1;

			points = static_cast<Point<dimensions>*>(std::realloc(points, maxPointCount * sizeof(Point<dimensions>)));
		}

		points[pointCount] = point;

		return pointCount++;
	}

	// Used as entry by the handle to remove a point.
	void removePoint(const uint32_t id, const float weighting) {
		// It possible tht the point for the handle is no longer stored here and has been moved downwards to some child on splitting. 
		// Thus we check for this and search the children if necesary. Note that the tree is guaranteed to be relatively shallow by the maxDepth parameter
		// meaning it should remain fast even in the worst cases. Additionally the majority of points will be stored at the bottom of the tree and never split.
		if (!pointHandle.isEmpty()) {

			pointHandle.removePoint(weighting);

			return; // Removal will be handled by upwards propagation.
		}

		removePoint(points[id], weighting); // Propagate removal upwards.

		freePointIndices.push_back(id);
	}

	// Used to propagate upwards point removal.
	void removePoint(const Point<dimensions>& point, const float weighting) {

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

		if (parent != nullptr) {
			parent->removePoint(point, weighting);
		}
	}

	// Calculates the repulsive energy omitting the constant multipliers.
	float calculateUnscaledRepulsiveEnergy(const Point<dimensions>& point, const float repulsiveExponent) {
		float squaredDistance = Point<dimensions>::squaredDistance(point, centreOfMass);

		// A point may not repel itself.
		if (squaredDistance == 0.0f) {
			return 0.0f;
		}

		// If we can consider this branch as a single mass. Leaf nodes may also be considered as a single mass. 
		// Leaf nodes may either be at maximum depth or contain exactly one point.
		if (squaredDistance >= maxWidthSquared || depth == maxDepth || maxPointCount == 1) {
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

		for (int i = 0; i < childCount; ++i) {
			if (childCreated[i]) {
				energy += children[i].calculateUnscaledRepulsiveEnergy(point, repulsiveExponent);
			}
		}

		return energy;
	}

	// Calculates the repulsive direction omitting the constant multipliers.
	Point<dimensions> calculateUnscaledRepulsiveDirection(const Point<dimensions>& point, const float repulsiveExponent) {
		float squaredDistance = Point<dimensions>::squaredDistance(point, centreOfMass);

		if (squaredDistance == 0.0f) {
			return Point<dimensions>();
		}

		// If we can consider this branch as a single mass. Leaf nodes may also be considered as a single mass. 
		// Leaf nodes may either be at maximum depth or contain exactly one point.
		if (squaredDistance >= maxWidthSquared || depth == maxDepth || maxPointCount == 1) {
			Point<dimensions> direction;

			for (int i = 0; i < dimensions; ++i) {
				direction[i] = (point[i] - centreOfMass[i]);
			}

			// We use the squared distance instead of the normal distance and thus add a factor 0.5f to the exponent.
			direction.scale(weight * std::powf(squaredDistance, 0.5f * repulsiveExponent - 1.0f));

			return direction;
		}

		// Otherwise this must not be a leaf node or be treade as a single mass and we must recurse into all child trees..
		Point<dimensions> direction = Point<dimensions>();

		for (int i = 0; i < childCount; ++i) {
			if (childCreated[i]) {
				direction.add(children[i].calculateUnscaledRepulsiveDirection(point, repulsiveExponent));
			}
		}

		return direction;
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
template <uint32_t dimensions, uint32_t maxDepth>
class QuadtreePointHandle {
public:
	void removePoint(const float weighting) {
		container->removePoint(id, weighting);
	}

private:
	Quadtree<dimensions, maxDepth>* container;
	uint32_t id;

	QuadtreePointHandle(Quadtree<dimensions, maxDepth>* container, uint32_t id) : container(container), id(id) {

	}

	// Empty handle constructor.
	QuadtreePointHandle() : container(nullptr), id(0) {

	}

	bool isEmpty() {
		return container == nullptr;
	}

	friend Quadtree<dimensions, maxDepth>;
};
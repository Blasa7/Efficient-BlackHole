#pragma once

#include <cstdint>
#include <cmath>

// Simple wrapper for points. At least, for 3 and 4 dimensional points the compiler will optimize the functions to use SIMD instructions where appropiate.
template <uint32_t dimensions>
class Point {
public:
	Point() {
		for (int i = 0; i < dimensions; ++i) {
			this->operator[](i) = 0.0f;
		}
	}

	Point(const float* position) {
		for (int i = 0; i < dimensions; ++i) {
			this->operator[](i) = position[i];
		}
	}

	const float* data() const {
		return position;
	}

	// Sum to current point.
	void add(const Point<dimensions>& point) {
		for (int i = 0; i < dimensions; ++i) {
			position[i] += point[i];
		}
	}

	// Scale current point.
	void scale(float scale) {
		for (int i = 0; i < dimensions; ++i) {
			position[i] *= scale;
		}
	}

	static inline float squaredDistance(const Point<dimensions>& a, const Point<dimensions>& b) {
		float sum = 0.0f;

		for (int i = 0; i < dimensions; ++i) {
			float diff = b[i] - a[i];
			sum += diff * diff;
		}

		return sum;
	}

	static inline float distance(const Point<dimensions>& a, const Point<dimensions>& b) {
		return std::sqrtf(squaredDistance(a, b));
	}

	inline float& operator[](size_t i) { return position[i]; }
	inline const float& operator[](std::size_t i) const { return position[i]; }

protected:
	// Position values, can be indexed by dimension.
	alignas(16) float position[dimensions];
};
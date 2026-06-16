#pragma once

#include <random>

// Helper for generating fairly good quality random numbers.
namespace random {
	thread_local std::mt19937 rng(std::random_device{}()); // Thread local to allow multithreaded use.

	inline int range(int min, int max) {
		return std::uniform_int_distribution<int>(min, max)(rng);
	}

	inline float nextf() {
		return std::uniform_real_distribution<float>()(rng);
	}

	inline double nextd() {
		return std::uniform_real_distribution<double>()(rng);
	}
}
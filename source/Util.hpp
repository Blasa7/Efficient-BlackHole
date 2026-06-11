#include <cmath>

#pragma once


//To calculate distance of two array
class CalculateUtil {
public :
	template <uint32_t dimensions>
	static double calcDist_DIM(double* a, double* b){
		double sum = 0.0;
		for (int i = 0; i < dimensions; i++){
			sum = sum + ((a[i] - b[i]) * (a[i] - b[i]));
		}
		return sqrt(sum);
	}
};

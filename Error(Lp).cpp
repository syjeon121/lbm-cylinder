#include "Error(Lp).h"
#include <iostream>
#include <cmath>
using namespace std;

float Error_Lp::Lp(int p, float *pArray1, float *pArray2, int n) {
	//p : the number of error
	//n : the number of node


	result = 0.0;
	sum = 0.0;
	int i;

	if (!p) cout << "Error! p is not netural number!" << endl;
		

	switch (p) {
	case 1:
		for (i = 0; i < n; i++) {
			sum = sum + abs(*pArray1 - *pArray2);

			pArray1 = pArray1 + 1;
			pArray2 = pArray2 + 1;
		}

		result = sum / n;
		break;

	case 2:
		for (i = 0; i < n; i++) {
			sum = sum + pow(abs(*pArray1 - *pArray2), 2);

			pArray1 = pArray1 + 1;
			pArray2 = pArray2 + 1;
		}

		result = sqrt(sum/n);
		break;

	case 33:
		float difference[3];
		float temp = 0.0;

		for (i = 0; i < n; i++) {
			difference[i] = abs(*pArray1 - *pArray2);
			if (difference[i] > temp) temp = difference[i];
		}
		
		result = temp;
		break;
	}

	
	
	return result;
}
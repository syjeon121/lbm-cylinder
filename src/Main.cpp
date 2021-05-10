#include "CPU.h"
#include "GPU.h"

int main() {

	int num;

	cout << "Choose the device!" << endl;
	cout << "1 : CPU" << endl;
	cout << "2 : GPU" << endl;
	cout << ">>";
	cin >> num;

	switch (num) {
	case 1:
		CPU(); break;
	case 2:
		GPU(); break;
	default:
		cout << "Error! Choose again"; return 0;
	}

	system("pause");
	return 0;
}
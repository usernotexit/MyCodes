#include<iostream>
#include<vector>
#include<Windows.h>
#include"Exercise.h"

int main() 
{
	char c;
	std::cout << "Welcome! Choose which subproject to show: (input: one number)\n"
		<< "1. Solve 84 linear equations by Gaussian elimination (1)\n"
		<< "2. Solve 100 * by * and advanced* (2.1)\n"
		<< "3. Solve 40'Hilbert equations by * and advanced* (2.2)\n"
		<< "4. Solve 100 * by Gaussian elimination (3.1)\n"
		<< "5. Solve 40'Hilbert equations by Gaussian elimination (3.2)\n"
		<< "6. I want them all!" << std::endl;
		
	std::cin >> c;
	cin.ignore(100, '\n');
	
	switch(c){
		case '1': exercise_1();break;
		case '2': exercise_2_1();break;
		case '3': exercise_2_2();break;
		case '4': exercise_3_1();break;
		case '5': exercise_3_2();break;
		case '6': {
			std::cout << "=================================================\n"
				<< "exercise 1" << std::endl;
			exercise_1();
			std::cout << "=================================================\n"
				<< "exercise 2.1" << std::endl;			
		 	exercise_2_1();
			std::cout << "=================================================\n"
				<< "exercise 2.2" << std::endl;
			exercise_2_2();
			std::cout << "=================================================\n"
				<< "exercise 3.1" << std::endl;
			exercise_3_1();
			std::cout << "=================================================\n"
				<< "exercise 3.2" << std::endl;
			exercise_3_2();
			break;
		}
	}
	
	return 0;
}

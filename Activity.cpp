
//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;
int main(int argc, char **argv){
	std::ifstream input("b657-wars.txt");
  	for(std::string line; getline(input, line);){
    		cout<<line<<endl;
    		}
    	}
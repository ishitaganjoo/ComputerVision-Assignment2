// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include <fstream>
#include <cmath>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

double euclideanDistance(SiftDescriptor image1_singleDesc, SiftDescriptor image2_singleDesc){
	double finalMagnitude = 0.0;
	for(int i=0; i<128; i++){
		 finalMagnitude += ((image1_singleDesc.descriptor[i]) - (image2_singleDesc.descriptor[i])) * ((image1_singleDesc.descriptor[i]) - (image2_singleDesc.descriptor[i]));
	}
	return sqrt(finalMagnitude);
}

int findNearestNeighbour(SiftDescriptor image1_singleDesc, vector<SiftDescriptor> imageDesc_2){
	double minDist = 1e6;
	double secondMin = 1e6;
	double threshold = 0.79;
	int neighborIndex = -1;
	for(int i=0; i<imageDesc_2.size(); i++){
		double distance = euclideanDistance(image1_singleDesc, imageDesc_2[i]);

		if(distance<minDist)
		{
			secondMin = minDist;
			minDist = distance;
			neighborIndex = i;
		}
		else if(distance<secondMin)
		{
			secondMin = distance;
		}

	}
	double ratio = minDist / secondMin;
	cout<<"ratio"<<ratio<<endl;
	cout<<"neighbor index"<<neighborIndex;
	if(ratio < threshold)
	{
		cout<<"threshold"<<threshold<<endl;
		cout<<"ratio inside"<<ratio<<endl;
		return neighborIndex;
	}
	return -1;
}

void findMatches(vector<SiftDescriptor> imageDesc_1, vector<SiftDescriptor> imageDesc_2, CImg<double>* finalImage, double imageSize, vector<SiftDescriptor> matchedVector){
	const unsigned char color[] = {255, 0, 0};
	int index = 0;
	for(int i=0; i<imageDesc_1.size(); i++){
		int indexNeighbor = findNearestNeighbour(imageDesc_1[i] , imageDesc_2);
		cout<<"indexNeighbour"<<indexNeighbor<<endl;
		if(indexNeighbor != -1)
		{
			cout<<"inside if"<<endl;
			cout<<imageDesc_1[i].col<<endl;
			cout<<matchedVector[0].col<<endl;
			matchedVector[index].col = imageDesc_1[i].col;
			matchedVector[index].row = imageDesc_1[i].row;
			for(int l=0; l<128; l++){
				matchedVector[index].descriptor[l] = imageDesc_1[i].descriptor[l];
			}
			index++;
			finalImage->draw_line(imageDesc_1[i].col , imageDesc_1[i].row, imageDesc_2[indexNeighbor].col+imageSize, imageDesc_2[indexNeighbor].row, color);
		}
	}
}

vector<SiftDescriptor> calculateDescriptors(CImg<double> image, string inputFileName){
	CImg<double> gray = image.get_RGBtoHSI().get_channel(2);
		vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);

	for(int i=0; i<descriptors.size(); i++)
		{
			cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
			for(int l=0; l<128; l++)
				cout << descriptors[i].descriptor[l] << "," ;
			cout << ")" << endl;

			for(int j=0; j<5; j++)
				for(int k=0; k<5; k++)
		if(j==2 || k==2)
			for(int p=0; p<3; p++)
				image(descriptors[i].col+k-1, descriptors[i].row+j-1, 0, p)=0;

		}
		string fileName = "sift"+inputFileName+".png";
		image.get_normalize(0, 255).save(fileName.c_str());
		return descriptors;
}


int main(int argc, char **argv)
{
  try {

    if(argc < 2)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    a2 part_id ..." << endl;
	return -1;
      }

    string part = argv[1];
    string inputFile = argv[2];
		string inputFile_2 = argv[3];

    if(part == "part1")
      {
				// This is just a bit of sample code to get you started, to
				// show how to use the SIFT library.

				CImg<double> input_image(inputFile.c_str());
				CImg<double> input_image_2(inputFile_2.c_str());

				CImg<double> finalImage = (input_image.get_append(input_image_2, 'x')).get_normalize(0, 255);

				vector<SiftDescriptor> desc1 = calculateDescriptors(input_image, "A");
				vector<SiftDescriptor> desc2 = calculateDescriptors(input_image_2, "B");

				//cout<<desc1[0]<<endl;

				vector<SiftDescriptor> matchedVector;

				findMatches(desc1, desc2, &finalImage, input_image.width(), matchedVector);

				finalImage.get_normalize(0, 255).save("final_2.png");

				// // convert image to grayscale
				// CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
				// 	vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
				//
				// for(int i=0; i<descriptors.size(); i++)
				//   {
				//     cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
				//     for(int l=0; l<128; l++)
				//       cout << descriptors[i].descriptor[l] << "," ;
				//     cout << ")" << endl;
				//
				//     for(int j=0; j<5; j++)
				//       for(int k=0; k<5; k++)
				// 	if(j==2 || k==2)
				// 	  for(int p=0; p<3; p++)
				// 	    input_image(descriptors[i].col+k-1, descriptors[i].row+j-1, 0, p)=0;
				//
				//   }
				//
				// 	CImg<double> input_image_2(inputFile_2.c_str());
				// 	CImg<double> gray_2 = input_image_2.get_RGBtoHSI().get_channel(2);
				//
				// 	vector<SiftDescriptor> descriptors_2 = Sift::compute_sift(gray_2);
				//
				// 	for(int i=0; i<descriptors_2.size(); i++)
				// 		{
				// 			cout << "Descriptor #" << i << ": x=" << descriptors_2[i].col << " y=" << descriptors[i]_2.row << " descriptor=(";
				// 			for(int l=0; l<128; l++)
				// 				cout << descriptors_2[i].descriptor[l] << "," ;
				// 			cout << ")" << endl;
				//
				// 			for(int j=0; j<5; j++)
				// 				for(int k=0; k<5; k++)
				// 		if(j==2 || k==2)
				// 			for(int p=0; p<3; p++)
				// 				input_image_2(descriptors_2[i].col+k-1, descriptors_2[i].row+j-1, 0, p)=0;
				//
				// 		}
				//
				//
				// 	input_image_2.get_normalize(0, 255).save("sift_2.png");
				//
				// input_image.get_normalize(0,255).save("sift.png");
      }
    else if(part == "part2")
      {
				CImg<double> input_image(inputFile.c_str());
				CImg<double> input_image_2(inputFile_2.c_str());

				vector<SiftDescriptor> desc1 = calculateDescriptors(input_image, "A");
				vector<SiftDescriptor> desc2 = calculateDescriptors(input_image_2, "B");

				CImg<double> finalImage = (input_image.get_append(input_image_2, 'x')).get_normalize(0, 255);

				vector<SiftDescriptor> matchedVector = vector<SiftDescriptor>();

				findMatches(desc1, desc2, &finalImage, input_image.width(), matchedVector);

				cout<<desc1.size()<<endl;

				int firstRandDesc = rand() % desc1.size();
				int secondRandDesc = rand() % desc1.size();
				int thirdRandDesc = rand() % desc1.size();
				int fourthRandDesc = rand() % desc1.size();

				cout<<firstRandDesc<<endl;
	// do something here!
      }
    else if(part == "part3")
      {
	// do something here!
      }
    else if(part == "part4"){
	cout<<"rohil"<<endl;
	std::ifstream input("b657-wars.txt");
  	for(std::string line; getline(input, line);){
    		cout<<line<<endl;
  	}
	}
    else
      throw std::string("unknown part!");

    // feel free to add more conditions for other parts (e.g. more specific)
    //  parts, for debugging, etc.
  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}

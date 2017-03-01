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
#include "randPoints.h"

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

void findMatches(vector<SiftDescriptor> imageDesc_1, vector<SiftDescriptor> imageDesc_2, CImg<double>* finalImage, double imageSize, vector<randPoints>* randomPoints){
	const unsigned char color[] = {255, 0, 0};
	for(int i=0; i<imageDesc_1.size(); i++){
		int indexNeighbor = findNearestNeighbour(imageDesc_1[i] , imageDesc_2);
		cout<<"indexNeighbour"<<indexNeighbor<<endl;
		if(indexNeighbor != -1)
		{
			randPoints obj = randPoints(imageDesc_1[i].col, imageDesc_1[i].row, imageDesc_2[indexNeighbor].col, imageDesc_2[indexNeighbor].row);
			//randPoints obj;
			//obj = randPoints(1,2,3,4);
			randomPoints->push_back(obj);
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

				vector<randPoints> randomPoints;

				findMatches(desc1, desc2, &finalImage, input_image.width(), &randomPoints);

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

				vector<randPoints> randomPoints;

				findMatches(desc1, desc2, &finalImage, input_image.width(), &randomPoints);

				for(int i=0; i<1; i++)
				{
					int firstRandDesc = rand() % randomPoints.size();
					int secondRandDesc = rand() % randomPoints.size();
					int thirdRandDesc = rand() % randomPoints.size();
					int fourthRandDesc = rand() % randomPoints.size();

					//vector<float> vectorRows;
					//vector<float> vectorCols;
					//vector<vector<float> > array_2d(8, vector<float>(9, 0));
					float x1 = randomPoints[firstRandDesc]._x;
					float y1 = randomPoints[firstRandDesc]._y;
					float x1P = randomPoints[firstRandDesc]._xP;
					float y1P = randomPoints[firstRandDesc]._yP;
				
					float x2 = randomPoints[secondRandDesc]._x;
					float y2 = randomPoints[secondRandDesc]._y;
					float x2P = randomPoints[secondRandDesc]._xP;
					float y2P = randomPoints[secondRandDesc]._yP;

					float x3 = randomPoints[thirdRandDesc]._x;
					float y3 = randomPoints[thirdRandDesc]._y;
					float x3P = randomPoints[thirdRandDesc]._xP;
					float y3P = randomPoints[thirdRandDesc]._yP;
				
					float x4 = randomPoints[fourthRandDesc]._x;
					float y4 = randomPoints[fourthRandDesc]._y;
					float x4P = randomPoints[fourthRandDesc]._xP;
					float y4P = randomPoints[fourthRandDesc]._yP;
					float array_2d[8][9] = {{-x1, -y1, -1, 0, 0, 0, x1*x1P, y1*x1P, x1P},
						    {0, 0, 0, -x1, -y1, -1, x1*y1P, y1*y1P, y1P},
						    {-x2, -y2, -1, 0, 0, 0, x2*x2P, y2*x2P, x2P},
						    {0, 0, 0, -x2, -y2, -1, x2*y2P, y2*y2P, y2P},
						    {-x3, -y3, -1, 0, 0, 0, x3*x3P, y3*x3P, x3P},
						    {0, 0, 0, -x3, -y3, -1, x3*y3P, y3*y3P, y3P},
						    {-x4, -y4, -1, 0, 0, 0, x4*x4P, y4*x4P, x4P},
						    {0, 0, 0, -x4, -y4, -1, x4*y4P, y4*y4P, y4P}};
					cout<<array_2d[0][0]<<"ishita"<<endl;
			}

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

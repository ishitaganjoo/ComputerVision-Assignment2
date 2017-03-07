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
#include <map>
#include <cstdlib>
#include <limits>
#include <random>
#include "randPoints.h"
#define kSize 28
#define wScalar 6

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
	double threshold = 0.75;
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
	//cout<<"ratio"<<ratio<<endl;
	//cout<<"neighbor index"<<neighborIndex;
	if(ratio < threshold)
	{
		//cout<<"threshold"<<threshold<<endl;
		//cout<<"ratio inside"<<ratio<<endl;
		return neighborIndex;
	}
	return -1;
}

int findMatches(vector<SiftDescriptor> imageDesc_1, vector<SiftDescriptor> imageDesc_2, CImg<double>* finalImage, double imageSize, vector<randPoints>* randomPoints,int countNoOfMatches, bool checkFlow){
	const unsigned char color[] = {255, 0, 0};
	for(int i=0; i<imageDesc_1.size(); i++){
		int indexNeighbor = findNearestNeighbour(imageDesc_1[i] , imageDesc_2);
		//cout<<"indexNeighbour"<<indexNeighbor<<endl;
		if(indexNeighbor != -1)
		{
			countNoOfMatches++;
			randPoints obj = randPoints(imageDesc_1[i].col, imageDesc_1[i].row, imageDesc_2[indexNeighbor].col, imageDesc_2[indexNeighbor].row);
			//randPoints obj;
			//obj = randPoints(1,2,3,4);
			randomPoints->push_back(obj);
			if(!checkFlow)
      		{
       
			finalImage->draw_line(imageDesc_1[i].col , imageDesc_1[i].row, imageDesc_2[indexNeighbor].col+imageSize, imageDesc_2[indexNeighbor].row, color);
			}
		}
	}
	return countNoOfMatches;
}

vector<SiftDescriptor> calculateDescriptors(CImg<double> image, string inputFileName){
	CImg<double> gray = image.get_RGBtoHSI().get_channel(2);
		vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);

	for(int i=0; i<descriptors.size(); i++)
		{
			//cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
			//for(int l=0; l<128; l++)
			//	cout << descriptors[i].descriptor[l] << "," ;
			//cout << ")" << endl;

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

//Question1 part2
int findBestMatchForImage(vector<SiftDescriptor> queryImage, vector<SiftDescriptor> inputImage, CImg<double>* finalImage, double imageSize, vector<randPoints>* randomPoints)
{
    int countNoOfMatches = 0;
    countNoOfMatches = findMatches(queryImage, inputImage, finalImage, imageSize, randomPoints, countNoOfMatches,true);
    //cout<<"inside bestMatches"<<countNoOfMatches<<endl;
    return countNoOfMatches;
}

vector<float> calculateHeuristicVector(SiftDescriptor desc, default_random_engine generator){

	vector<float> fV(kSize);
	vector<float> x(128);

	normal_distribution<double> distribution(0.0, 1.0);

    for(int j = 0; j < 128; j++){
        x[j] = distribution(generator);
    }
	for(int i=0; i<kSize; i++){
		for(int k=0; k<128; k++){
			fV[i] += desc.descriptor[k]*x[k];
			fV[i] = floor(fV[i] / wScalar);
		}
	}

	return fV;
}

void matchHeuristicVectors(vector<float> heuristicImage1, SiftDescriptor image1128D, vector< vector<float> > heuristicAllImage2, vector<SiftDescriptor> image2128D, CImg<double>* finalImage, double imageSize){
	int newCount = 0;
	double minDist = 1e6;
	double secondMin = 1e6;
	double threshold = 0.75;
	int neighborIndex = -1;
	for(int i=0; i<heuristicAllImage2.size(); i++){
		double distance = 100000;
		int count = 0;
		for(int k=0; k<kSize; k++){
			if(heuristicImage1[k] == heuristicAllImage2[i][k]){
				count++;
			}
			else{
				count--;
			}
		}
		if(count == kSize){
			distance = euclideanDistance(image1128D, image2128D[i]);

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
	}
	double ratio = minDist / secondMin;
	if(ratio < threshold)
	{
		if(neighborIndex != -1)
		{

			//  randPoints obj = randPoints(imageDesc_1[i].col, imageDesc_1[i].row, imageDesc_2[indexNeighbor].col, imageDesc_2[indexNeighbor].row);
			// //randPoints obj;
			// //obj = randPoints(1,2,3,4);
			// randomPoints->push_back(obj);
			const unsigned char color[] = {255, 0, 0};
			finalImage->draw_line(image1128D.col , image1128D.row, image2128D[neighborIndex].col+imageSize, image2128D[neighborIndex].row, color);
		}
	}
}


// Calculate min homography matrix.
CImg<float> calculateMinHomographyMatrix(vector<SiftDescriptor> desc1, vector<SiftDescriptor> desc2, int inputWidth)
{		
		
        CImg<double> finalImage;
		vector<randPoints> randomPoints;

		findMatches(desc1, desc2, &finalImage, inputWidth, &randomPoints, 0, true);
        
		CImg<float> points(2, 4);
		double minErrorDist = 100000000.0;
		CImg<float> minHomographyMatrix(3,3);
		for(int i = 0 ; i < 3; i++)
			{
			for(int j = 0; j < 3; j++)
				{
					minHomographyMatrix(i,j) = 0.0;
				}
			}
			
			for(int i=0 ; i<2000 ; i++)
				{		
					int firstRandDesc = rand() % randomPoints.size();
					int secondRandDesc = rand() % randomPoints.size();
					int thirdRandDesc = rand() % randomPoints.size();
					int fourthRandDesc = rand() % randomPoints.size();

					CImg<float> origTranslatedPoints(2, 4);
					
					float x1 = randomPoints[firstRandDesc]._x;
					points(0, 0) = randomPoints[firstRandDesc]._x;
					float y1 = randomPoints[firstRandDesc]._y;
					points(1, 0) = randomPoints[firstRandDesc]._y;
					float x1P = randomPoints[firstRandDesc]._xP;
					origTranslatedPoints(0, 0) = randomPoints[firstRandDesc]._xP;
					float y1P = randomPoints[firstRandDesc]._yP;
					origTranslatedPoints(1, 0) = randomPoints[firstRandDesc]._yP;
						
					float x2 = randomPoints[secondRandDesc]._x;
					points(0, 1) = randomPoints[secondRandDesc]._x;
					float y2 = randomPoints[secondRandDesc]._y;
					points(1, 1) = randomPoints[secondRandDesc]._y;
					float x2P = randomPoints[secondRandDesc]._xP;
					origTranslatedPoints(0, 1) = randomPoints[secondRandDesc]._xP;
					float y2P = randomPoints[secondRandDesc]._yP;
					origTranslatedPoints(1, 1) = randomPoints[secondRandDesc]._yP;
					
					float x3 = randomPoints[thirdRandDesc]._x;
					points(0, 2) = randomPoints[thirdRandDesc]._x;
					float y3 = randomPoints[thirdRandDesc]._y;
					points(1, 2) = randomPoints[thirdRandDesc]._y;
					float x3P = randomPoints[thirdRandDesc]._xP;
					origTranslatedPoints(0, 2) = randomPoints[thirdRandDesc]._xP;
					float y3P = randomPoints[thirdRandDesc]._yP;
					origTranslatedPoints(1, 2) = randomPoints[thirdRandDesc]._yP;
				
					float x4 = randomPoints[fourthRandDesc]._x;
					points(0, 3) = randomPoints[fourthRandDesc]._x;
					float y4 = randomPoints[fourthRandDesc]._y;
					points(1, 3) = randomPoints[fourthRandDesc]._y;
					float x4P = randomPoints[fourthRandDesc]._xP;
					origTranslatedPoints(0, 3) = randomPoints[fourthRandDesc]._xP;
					float y4P = randomPoints[fourthRandDesc]._yP;
					origTranslatedPoints(1, 3) = randomPoints[fourthRandDesc]._yP;	
					CImg<float> array_2d(8,8);
					for(int i=0; i<8; i++){
						for(int j=0; j<8; j++){
							if(i % 2 == 0){
								if(j == 3 || j == 4 || j == 5){
									array_2d(i, j) = 0;
								}
								if(j == 2){
									array_2d(i, j) = 1;
								}
							}
							if(i % 2 != 0){
								if(j == 0 || j == 1 || j == 2){
									array_2d(i, j) = 0;
								}
								if(j == 5){
									array_2d(i, j) = 1;
								}
							}
						}
					}
					
					array_2d(0,0) = x1;
					array_2d(0, 1) = y1;
					array_2d(0, 6) = -(x1*x1P);
					array_2d(0, 7) = -(y1*x1P);
					array_2d(1, 3) = x1;
					array_2d(1, 4) = y1;
					array_2d(1, 6) = -(x1*y1P);
					array_2d(1, 7) = -(y1*y1P);

					array_2d(2,0) = x2;
					array_2d(2, 1) = y2;
					array_2d(2, 6) = -(x2*x2P);
					array_2d(2, 7) = -(y2*x2P);
					array_2d(3, 3) = x2;
					array_2d(3, 4) = y2;
					array_2d(3, 6) = -(x2*y2P);
					array_2d(3, 7) = -(y2*y2P);

					array_2d(4,0) = x3;
					array_2d(4, 1) = y3;
					array_2d(4, 6) = -(x3*x3P);
					array_2d(4, 7) = -(y3*x3P);
					array_2d(5, 3) = x3;
					array_2d(5, 4) = y3;
					array_2d(5, 6) = -(x3*y3P);
					array_2d(5, 7) = -(y3*y3P);

					array_2d(6,0) = x4;
					array_2d(6, 1) = y4;
					array_2d(6, 6) = -(x4*x4P);
					array_2d(6, 7) = -(y4*x4P);
					array_2d(7, 3) = x4;
					array_2d(7, 4) = y4;
					array_2d(7, 6) = -(x4*y4P);
					array_2d(7, 7) = -(y4*y4P);

					//float d = determinantOfMatrix(array_2d, 8);
					CImg<float> inverseMat(8, 8);
					//inverse(array_2d, inverseMat);
					//cout<<"inverse mat"<<inverseMat(0,0)<<endl;
					inverseMat = array_2d.invert();
					//cout<<inverseMat(0,0)<<"invert"<<endl;

					vector<float> transCoord(8);
					transCoord[0] = x1P;
					transCoord[1] = y1P;
					transCoord[2] = x2P;
					transCoord[3] = y2P;
					transCoord[4] = x3P;
					transCoord[5] = y3P;
					transCoord[6] = x4P;
					transCoord[7] = y4P;

					vector<float> homography(9);
					for(int i = 0 ; i < 9; i++)
					{
						homography[i] = 0.0;
					}

					for(int i=0; i<inverseMat.height(); i++){
						for(int j=0; j<inverseMat.width(); j++){
							homography[i] += inverseMat(i, j) * transCoord[j];
						}
					}
					homography[8] = 1;
					
					CImg<float> homographyMatrix(3,3);
					for(int i = 0 ; i < 3; i++)
					{
						for(int j = 0; j < 3; j++)
						{
							homographyMatrix(i,j) = 0.0;
						}

					}
					int index = 0;
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							homographyMatrix(i, j) = homography[index];
							index++;
						}
					}

					CImg<float> translatedCoord(2, 4);

					int pointIndex = 0;

					while(pointIndex < 4){
						translatedCoord(0, pointIndex) = homographyMatrix(0, 0)*points(0,pointIndex) + homographyMatrix(0, 1)*points(1, pointIndex) + homographyMatrix(0, 2);
						translatedCoord(1, pointIndex) = homographyMatrix(1, 0)*points(0, pointIndex) + homographyMatrix(1, 1)*points(1, pointIndex) + homographyMatrix(1, 2);
						pointIndex++;
					}

					double errorDistance = 0.0;

					pointIndex = 0;
					while(pointIndex < 4){
						errorDistance += sqrt(pow(translatedCoord(0, pointIndex) - origTranslatedPoints(0, pointIndex),2) + pow(translatedCoord(1,pointIndex) - origTranslatedPoints(1, pointIndex),2));
						pointIndex++;
					}

					if(errorDistance < minErrorDist)
					{
						minErrorDist = errorDistance;
						minHomographyMatrix = homographyMatrix ;
					}
				}
				
				return minHomographyMatrix;
}

// Warp my image
void warpImage(CImg<double> inverted, CImg<double> input_image, string filename)
{
	    int height = input_image.height();
		int width = input_image.width();
		CImg<double> output_image(width, height, 1, 3, 0.0); 
        for (int i = 0; i < width; i++)
				{  
				for (int j = 0; j < height; j++)
					{		
						double x = inverted(0, 0) * i + inverted(0, 1) * j + inverted(0, 2);
						double y = inverted(1, 0) * i + inverted(1, 1) * j + inverted(1, 2);
						double z = inverted(2, 0) * i + inverted(2, 1) * j + inverted(2, 2);
						if (z>0){
							x /= z;
							y /= z;
			  				   }
						if (x>=0 && x<width && y>=0 && y<height) 
	  						{
							output_image(i, j, 0) = input_image(x, y, 0);
							output_image(i, j, 1) = input_image(x, y, 1);
							output_image(i, j, 2) = input_image(x, y, 2);
							}
						else
							{
							output_image(i, j, 0) = 0.0;
							output_image(i, j, 1) = 0.0;
							output_image(i, j, 2) = 0.0;
							}
					}
			    }
		string name = filename+"-warped.png";		    
	    output_image.save(name.c_str());	

}

CImg<double> getInverseTransformMatrix()
{   CImg<double> inverseTransform(3, 3);
		inverseTransform(0, 0) = 0.907;
		inverseTransform(0, 1) = 0.258;
		inverseTransform(0, 2) =  -182;
		inverseTransform(1, 0) =  -0.153;
		inverseTransform(1, 1) =  1.44;
		inverseTransform(1, 2) =  58;
		inverseTransform(2, 0) = -0.000306;
		inverseTransform(2, 1) = 0.000731;
		inverseTransform(2, 2) = 1;
		return inverseTransform;
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

    if(part == "part1" && argc == 4)
      {
				// This is just a bit of sample code to get you started, to
				// show how to use the SIFT library.
      			string inputFile_2 = argv[3];
				CImg<double> input_image(inputFile.c_str());
				CImg<double> input_image_2(inputFile_2.c_str());

				CImg<double> finalImage = (input_image.get_append(input_image_2, 'x')).get_normalize(0, 255);

				vector<SiftDescriptor> desc1 = calculateDescriptors(input_image, "A");
				vector<SiftDescriptor> desc2 = calculateDescriptors(input_image_2, "B");

				//cout<<desc1[0]<<endl;

				vector<randPoints> randomPoints;

				findMatches(desc1, desc2, &finalImage, input_image.width(), &randomPoints, 0, false);

				finalImage.get_normalize(0, 255).save("final_2.png");

      }
			else if(part == "part1" && argc > 4)
			{

					//create a map
					map<int,string> countMap;
					map<int,string>::const_iterator it;
					map<int,string>::const_iterator it1;
					int length = argc-3 ; //length of args
					string inputFile_1 = argv[2];
					CImg<double> query_image(inputFile_1.c_str());
					//cout<<length<<"length is"<<endl;
					vector<SiftDescriptor> desc1 = calculateDescriptors(query_image, "A");

					//start the for loop
					for(int i = 0; i < length; i++)
					{
					string inputFile_match = argv[i+3];
					CImg<double> input_image_2(inputFile_match.c_str());

					CImg<double> finalImage = (query_image.get_append(input_image_2, 'x')).get_normalize(0, 255);


					vector<SiftDescriptor> desc2 = calculateDescriptors(input_image_2, "B");

					vector<randPoints> randomPoints;
					int count=0;
					count = findBestMatchForImage(desc1, desc2, &finalImage, query_image.width(), &randomPoints);

					countMap[count] =  argv[i+3]; //name of file;
                    }

					it = countMap.end();
					it--;
					it1 = countMap.begin();
					it1--;
					for(it; it!=it1; it--)
					{
							//cout<<"count order is"<<it->first<<endl;
							cout<<"image order is"<<it->second<<endl;
					}

			}
    else if(part == "part2")
      {
      			string inputFile_2 = argv[3];
				CImg<double> input_image(inputFile.c_str());
				CImg<double> input_image_2(inputFile_2.c_str());

				double size = input_image.width();

				vector<SiftDescriptor> desc1 = calculateDescriptors(input_image, "A");
				vector<SiftDescriptor> desc2 = calculateDescriptors(input_image_2, "B");

				CImg<double> finalImage = (input_image.get_append(input_image_2, 'x')).get_normalize(0, 255);

				vector<randPoints> randomPoints;

				findMatches(desc1, desc2, &finalImage, input_image.width(), &randomPoints, 0, false);

				//map<double, vector<float> > errorCal;
				CImg<float> points(2, 4);

				double minErrorDist = 100000000.0;
				CImg<float> minHomographyMatrix(3,3);
				for(int i = 0 ; i < 3; i++)
				{
					for(int j = 0; j < 3; j++)
					{
						minHomographyMatrix(i,j) = 0.0;
					}

				}
				for(int i=0 ; i<2000 ; i++)
				{
					int firstRandDesc = rand() % randomPoints.size();
					int secondRandDesc = rand() % randomPoints.size();
					int thirdRandDesc = rand() % randomPoints.size();
					int fourthRandDesc = rand() % randomPoints.size();

					CImg<float> origTranslatedPoints(2, 4);

					float x1 = randomPoints[firstRandDesc]._x;
					points(0, 0) = randomPoints[firstRandDesc]._x;
					float y1 = randomPoints[firstRandDesc]._y;
					points(1, 0) = randomPoints[firstRandDesc]._y;
					float x1P = randomPoints[firstRandDesc]._xP;
					origTranslatedPoints(0, 0) = randomPoints[firstRandDesc]._xP;
					float y1P = randomPoints[firstRandDesc]._yP;
					origTranslatedPoints(1, 0) = randomPoints[firstRandDesc]._yP;

					float x2 = randomPoints[secondRandDesc]._x;
					points(0, 1) = randomPoints[secondRandDesc]._x;
					float y2 = randomPoints[secondRandDesc]._y;
					points(1, 1) = randomPoints[secondRandDesc]._y;
					float x2P = randomPoints[secondRandDesc]._xP;
					origTranslatedPoints(0, 1) = randomPoints[secondRandDesc]._xP;
					float y2P = randomPoints[secondRandDesc]._yP;
					origTranslatedPoints(1, 1) = randomPoints[secondRandDesc]._yP;

					float x3 = randomPoints[thirdRandDesc]._x;
					points(0, 2) = randomPoints[thirdRandDesc]._x;
					float y3 = randomPoints[thirdRandDesc]._y;
					points(1, 2) = randomPoints[thirdRandDesc]._y;
					float x3P = randomPoints[thirdRandDesc]._xP;
					origTranslatedPoints(0, 2) = randomPoints[thirdRandDesc]._xP;
					float y3P = randomPoints[thirdRandDesc]._yP;
					origTranslatedPoints(1, 2) = randomPoints[thirdRandDesc]._yP;

					float x4 = randomPoints[fourthRandDesc]._x;
					points(0, 3) = randomPoints[fourthRandDesc]._x;
					float y4 = randomPoints[fourthRandDesc]._y;
					points(1, 3) = randomPoints[fourthRandDesc]._y;
					float x4P = randomPoints[fourthRandDesc]._xP;
					origTranslatedPoints(0, 3) = randomPoints[fourthRandDesc]._xP;
					float y4P = randomPoints[fourthRandDesc]._yP;
					origTranslatedPoints(1, 3) = randomPoints[fourthRandDesc]._yP;

					CImg<float> array_2d(8,8);
					for(int i=0; i<8; i++){
						for(int j=0; j<8; j++){
							if(i % 2 == 0){
								if(j == 3 || j == 4 || j == 5){
									array_2d(i, j) = 0;
								}
								if(j == 2){
									array_2d(i, j) = 1;
								}
							}
							if(i % 2 != 0){
								if(j == 0 || j == 1 || j == 2){
									array_2d(i, j) = 0;
								}
								if(j == 5){
									array_2d(i, j) = 1;
								}
							}
						}
					}
					array_2d(0,0) = x1;
					array_2d(0, 1) = y1;
					array_2d(0, 6) = -(x1*x1P);
					array_2d(0, 7) = -(y1*x1P);
					array_2d(1, 3) = x1;
					array_2d(1, 4) = y1;
					array_2d(1, 6) = -(x1*y1P);
					array_2d(1, 7) = -(y1*y1P);

					array_2d(2,0) = x2;
					array_2d(2, 1) = y2;
					array_2d(2, 6) = -(x2*x2P);
					array_2d(2, 7) = -(y2*x2P);
					array_2d(3, 3) = x2;
					array_2d(3, 4) = y2;
					array_2d(3, 6) = -(x2*y2P);
					array_2d(3, 7) = -(y2*y2P);

					array_2d(4,0) = x3;
					array_2d(4, 1) = y3;
					array_2d(4, 6) = -(x3*x3P);
					array_2d(4, 7) = -(y3*x3P);
					array_2d(5, 3) = x3;
					array_2d(5, 4) = y3;
					array_2d(5, 6) = -(x3*y3P);
					array_2d(5, 7) = -(y3*y3P);

					array_2d(6,0) = x4;
					array_2d(6, 1) = y4;
					array_2d(6, 6) = -(x4*x4P);
					array_2d(6, 7) = -(y4*x4P);
					array_2d(7, 3) = x4;
					array_2d(7, 4) = y4;
					array_2d(7, 6) = -(x4*y4P);
					array_2d(7, 7) = -(y4*y4P);

					//float d = determinantOfMatrix(array_2d, 8);
					CImg<float> inverseMat(8, 8);
					//inverse(array_2d, inverseMat);
					//cout<<"inverse mat"<<inverseMat(0,0)<<endl;
					inverseMat = array_2d.invert();
					cout<<inverseMat(0,0)<<"invert"<<endl;

					vector<float> transCoord(8);
					transCoord[0] = x1P;
					transCoord[1] = y1P;
					transCoord[2] = x2P;
					transCoord[3] = y2P;
					transCoord[4] = x3P;
					transCoord[5] = y3P;
					transCoord[6] = x4P;
					transCoord[7] = y4P;

					vector<float> homography(9);
					for(int i = 0 ; i < 9; i++)
					{
						homography[i] = 0.0;
					}

					for(int i=0; i<inverseMat.height(); i++){
						for(int j=0; j<inverseMat.width(); j++){
							homography[i] += inverseMat(i, j) * transCoord[j];
						}
					}
					homography[8] = 1;

					CImg<float> homographyMatrix(3,3);
					for(int i = 0 ; i < 3; i++)
					{
						for(int j = 0; j < 3; j++)
						{
							homographyMatrix(i,j) = 0.0;
						}

					}
					int index = 0;
					for(int i=0; i<3; i++){
						for(int j=0; j<3; j++){
							homographyMatrix(i, j) = homography[index];
							index++;
						}
					}
					CImg<float> translatedCoord(2, 4);

					int pointIndex = 0;

					while(pointIndex < 4){
						translatedCoord(0, pointIndex) = homographyMatrix(0, 0)*points(0,pointIndex) + homographyMatrix(0, 1)*points(1, pointIndex) + homographyMatrix(0, 2);
						translatedCoord(1, pointIndex) = homographyMatrix(1, 0)*points(0, pointIndex) + homographyMatrix(1, 1)*points(1, pointIndex) + homographyMatrix(1, 2);
						pointIndex++;
					}

					double errorDistance = 0.0;

					pointIndex = 0;
					while(pointIndex < 4){
						errorDistance += sqrt(pow(translatedCoord(0, pointIndex) - origTranslatedPoints(0, pointIndex),2) + pow(translatedCoord(1,pointIndex) - origTranslatedPoints(1, pointIndex),2));
						pointIndex++;
					}
					cout<<"errorDist"<<errorDistance<<endl;
					//errorCal[errorDistance] = homography;
					if(errorDistance < minErrorDist)
					{
						minErrorDist = errorDistance;
						minHomographyMatrix = homographyMatrix ;
					}
			}

			cout<<minErrorDist<<" min error"<<endl;

			for(int i=0; i<3; i++){
				for(int j=0; j<3; j++){
					cout<<minHomographyMatrix(i, j)<<endl;
				}

			CImg<float> finalTranslatedCoord(2, randomPoints.size());

			int pointIndex = 0;

			const unsigned char color[] = {255, 0, 0};


			while(pointIndex < randomPoints.size()){
				finalTranslatedCoord(0, pointIndex) = minHomographyMatrix(0, 0)*randomPoints[pointIndex]._x + minHomographyMatrix(0, 1)*randomPoints[pointIndex]._y + minHomographyMatrix(0, 2);
				finalTranslatedCoord(1, pointIndex) = minHomographyMatrix(1, 0)*randomPoints[pointIndex]._x + minHomographyMatrix(1, 1)*randomPoints[pointIndex]._y + minHomographyMatrix(1, 2);
				finalImage.draw_line(randomPoints[pointIndex]._x , randomPoints[pointIndex]._y, finalTranslatedCoord(0, pointIndex)+size, finalTranslatedCoord(1, pointIndex), color);
				pointIndex++;
			}

			finalImage.get_normalize(0, 255).save("final_part2.png");


      }
}
			else if(part == "rohil"){
				string inputFile_2 = argv[3];
				CImg<double> input_image(inputFile.c_str());
				CImg<double> input_image_2(inputFile_2.c_str());

				CImg<double> finalImage = (input_image.get_append(input_image_2, 'x')).get_normalize(0, 255);

				vector<SiftDescriptor> desc1 = calculateDescriptors(input_image, "A");
				vector<SiftDescriptor> desc2 = calculateDescriptors(input_image_2, "B");

				vector< vector<float> > heuristicVecImage1(desc1.size());
				vector< vector<float> > heuristicVecImage2(desc2.size());

				for(int i=0; i<desc1.size(); i++){
					default_random_engine generator(time(0));
					heuristicVecImage1[i] = calculateHeuristicVector(desc1[i], generator);
					//matchHeuristicVectors(heuristicVec, desc1[i], desc2);
				}
				for(int i=0; i<desc2.size(); i++){
					default_random_engine generator(time(0));
					heuristicVecImage2[i] = calculateHeuristicVector(desc2[i], generator);
				}

				for(int i=0; i<heuristicVecImage1.size(); i++){
					matchHeuristicVectors(heuristicVecImage1[i], desc1[i], heuristicVecImage2, desc2, &finalImage, input_image.width());
				}
				finalImage.get_normalize(0, 255).save("final_part2_2.png");



			}
        else if(part == "part3")
      { // Question3/part1
      	if (argc == 3){
	  		CImg<float> input_image(inputFile.c_str());
			CImg<double> inverseTransform(3, 3); 
			inverseTransform = getInverseTransformMatrix();
			CImg<double> inverted(3, 3); 
			inverted = inverseTransform.invert();
			warpImage(inverted, input_image, inputFile);
	    }

	    // Question3/part2
	    else  {
	    	       //Iterate over the arguments
	    			int length = argc-3 ; //length of args
					string inputFile_1 = argv[2];
					CImg<double> query_image(inputFile_1.c_str());
					vector<SiftDescriptor> desc1 = calculateDescriptors(query_image, "A");
					int inputWidth = query_image.width();
					//start the for loop
					for(int i = 0; i < length; i++)
					{
					string inputFileWarp = argv[i+3];
					CImg<double> input_image_2(inputFileWarp.c_str());
					vector<SiftDescriptor> desc2 = calculateDescriptors(input_image_2, "B");
					//Calculate minimum homography matrix
					CImg<float> minHomographyMatrix(3,3);
					minHomographyMatrix = calculateMinHomographyMatrix(desc1, desc2, inputWidth);

					CImg<float> minHomographyMatrixInvert(3,3);
					minHomographyMatrixInvert = minHomographyMatrix.invert();
					int index = inputFileWarp.find_last_of("/"); 
					string filename = inputFileWarp.substr(index+1);
					int indexExt = filename.find_last_of(".");
					
					warpImage(minHomographyMatrixInvert, input_image_2, filename.substr(0, indexExt));
				}
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

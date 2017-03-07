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

void getCofactor(CImg<float> mat, CImg<float> temp, int p, int q, int n)
{
		int i = 0, j = 0;

		// Looping for each element of the matrix
		for (int row = 0; row < n; row++)
		{
				for (int col = 0; col < n; col++)
				{
						//  Copying into temporary matrix only those element
						//  which are not in given row and column
						if (row != p && col != q)
						{
								//temp[i][j++] = mat[row][col];
								temp(i, j++) = mat(row, col);

								// Row is filled, so increase row index and
								// reset col index
								if (j == n - 1)
								{
										j = 0;
										i++;
								}
						}
				}
		}
}

int determinantOfMatrix(CImg<float> mat, int n)
{
		int D = 0; // Initialize result

		// //  Base case : if matrix contains single element
		 if (n == 1)
		     return mat(0, 0);

		//int temp[N][N]; // To store cofactors
		CImg<float> temp(8,8);

		int sign = 1;  // To store sign multiplier

		 // Iterate for each element of first row
		for (int f = 0; f < n; f++)
		{
				// Getting Cofactor of mat[0][f]
				getCofactor(mat, temp, 0, f, n);
				D += sign * mat(0, f) * determinantOfMatrix(temp, n - 1);

				// terms are to be added with alternate sign
				sign = -sign;
		}

		return D;
}

void adjoint(CImg<float> A, CImg<float> adj)
{
    // if (N == 1)
    // {
    //     adj(0, 0) = 1;
    //     return;
    // }

    // temp is used to store cofactors of A[][]
    int sign = 1;
		CImg<float> temp;

    for (int i=0; i<8; i++)
    {
        for (int j=0; j<8; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, 8);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj(j, i) = (sign)*(determinantOfMatrix(temp, 7));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(CImg<float> A, CImg<float> inverseMat)
{
    // Find determinant of A[][]
    int det = determinantOfMatrix(A, 8);
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }

    // Find adjoint
    CImg<float> adj(8, 8);
    adjoint(A, adj);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<8; i++)
        for (int j=0; j<8; j++)
            inverseMat(i, j) = adj(i, j)/float(det);

    return true;
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

					CImg<float> points(2, 4);
					CImg<float> origTranslatedPoints(2, 4);

					// float x1 = randomPoints[firstRandDesc]._x;
					// float y1 = randomPoints[firstRandDesc]._y;
					// float x1P = randomPoints[firstRandDesc]._xP;
					// float y1P = randomPoints[firstRandDesc]._yP;
					//
					// float x2 = randomPoints[secondRandDesc]._x;
					// float y2 = randomPoints[secondRandDesc]._y;
					// float x2P = randomPoints[secondRandDesc]._xP;
					// float y2P = randomPoints[secondRandDesc]._yP;
					//
					// float x3 = randomPoints[thirdRandDesc]._x;
					// float y3 = randomPoints[thirdRandDesc]._y;
					// float x3P = randomPoints[thirdRandDesc]._xP;
					// float y3P = randomPoints[thirdRandDesc]._yP;
					//
					// float x4 = randomPoints[fourthRandDesc]._x;
					// float y4 = randomPoints[fourthRandDesc]._y;
					// float x4P = randomPoints[fourthRandDesc]._xP;
					// float y4P = randomPoints[fourthRandDesc]._yP;

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

						//	{{-x1, -y1, -1, 0, 0, 0, x1*x1P, y1*x1P, x1P},
						  //  {0, 0, 0, -x1, -y1, -1, x1*y1P, y1*y1P, y1P},
						   // {-x2, -y2, -1, 0, 0, 0, x2*x2P, y2*x2P, x2P},
						   // {0, 0, 0, -x2, -y2, -1, x2*y2P, y2*y2P, y2P},
						   // {-x3, -y3, -1, 0, 0, 0, x3*x3P, y3*x3P, x3P},
						   // {0, 0, 0, -x3, -y3, -1, x3*y3P, y3*y3P, y3P},
						   // {-x4, -y4, -1, 0, 0, 0, x4*x4P, y4*x4P, x4P},
						   // {0, 0, 0, -x4, -y4, -1, x4*y4P, y4*y4P, y4P}};
					//CImg<float> array_2d(8,9);
					//array_2d = {[-x1, -y1, -1, 0, 0, 0, x1*x1P, y1*x1P, x1P],
					//	    [0, 0, 0, -x1, -y1, -1, x1*y1P, y1*y1P, y1P],
					//	    [-x2, -y2, -1, 0, 0, 0, x2*x2P, y2*x2P, x2P],
					//	    [0, 0, 0, -x2, -y2, -1, x2*y2P, y2*y2P, y2P],
					//	    [-x3, -y3, -1, 0, 0, 0, x3*x3P, y3*x3P, x3P],
					//	    [0, 0, 0, -x3, -y3, -1, x3*y3P, y3*y3P, y3P],
					//	    [-x4, -y4, -1, 0, 0, 0, x4*x4P, y4*x4P, x4P],
					//	    [0, 0, 0, -x4, -y4, -1, x4*y4P, y4*y4P, y4P]};

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
					for(int i=0; i<inverseMat.height(); i++){
						for(int j=0; j<inverseMat.width(); j++){
							homography[i] += inverseMat(i, j) * transCoord[j];
						}
					}
					homography[8] = 1;
					for(int i=0; i<homography.size(); i++){
						cout<<"rohil "<<homography[i]<<endl;
					}
					CImg<float> homographyMatrix(3,3);
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

					// float imageCoordX = homographyMatrix(0,0)*x1 + homographyMatrix(0,1)*y1 + homographyMatrix(0, 2);
					// float imageCoordY = homographyMatrix(1,0)*x1 + homographyMatrix(1,1)*y1 + homographyMatrix(1, 2);
					cout<<"original "<<x1P<<endl;
					cout<<"tr "<<translatedCoord(0, 0)<<endl;
					cout<<"original "<<y1P<<endl;
					cout<<"tr "<<translatedCoord(1, 0)<<endl;
					cout<<"original "<<x2P<<endl;
					cout<<"tr "<<translatedCoord(0, 1)<<endl;
					cout<<"original "<<y2P<<endl;
					cout<<"tr "<<translatedCoord(1, 1)<<endl;

					double errorDistance = 0.0;
					pointIndex = 0;
					while(pointIndex < 4){
						errorDistance += sqrt(pow(translatedCoord(0, pointIndex) - origTranslatedPoints(0, pointIndex),2) - pow(translatedCoord(1,pointIndex) - origTranslatedPoints(1, pointIndex),2));
						pointIndex++;
					}
					cout<<"errorDist"<<errorDistance<<endl;


			}
      }
        else if(part == "part3")
      { // Question3/part1
      	if (argc == 2){
	  		CImg<float> input_image(inputFile.c_str());
			int height = input_image.height();
			int width = input_image.width();
			CImg<double> output_image(width, height, 1, 3, 0.0); 
			CImg<double> inverseTransform(3, 3); 
			inverseTransform = getInverseTransformMatrix();
			CImg<double> inverted(3, 3); 
			inverted = inverseTransform.invert();
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
	    output_image.save("warped_image.png");	
	    }
	    // Question3/part2
	    else {

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

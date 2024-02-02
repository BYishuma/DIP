#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include<string.h>
#include <ctype.h>
#include "proto.h"
#include <math.h>
#include<stdbool.h>
#include <vector>
#include <stack>
#define PI 3.141592

void SavePNMImage(Image*, char*);
Image* OTSU(Image* image);
double Findsigma(int k, double mg, double pi[]);
void Findpi(int image->Height, int image->Width, long histogram[], double pi[]);
double p1(int k, double pi[]);
double m(int k, double pi[]);
double abs(double a, double b);
void Findhistogram(int height, int width, Image* image, long histogram[]);
void Findpi(int height, int width, long histogram[], double pi[]);
Image* AverageFilter(Image* image, int size);
Image* MovingAverageThresholding(Image* image);
Image* PartialOtus(Image* image);
Image* RegionGrowing(Image* image, float threshold);

int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//septagon_noisy_shaded.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//spot_shaded_text_image.pgm ";
    char* input3 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//large_septagon_gaussian_noise_mean_0_std_50_added.pgm";
    char* input4 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//noisy_region.pgm";
    char* input5 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//defective_weld.pgm";
    char* output;

    //Task1
    //output = "output//large_septagon_gaussian_noise_mean_0_std_50_addedOtusAvg.pgm";
    //TestReadImage(input3, output);

    //Task2
    //output = "output//septagon_noisy_shadedPar.pgm";
    //TestReadImage(input1, output);

    //Task3
    output = "output//spot_shaded_text_imageAMT.pgm";
    TestReadImage(input2, output);
    
    //Task4
    //output= "output//defective_weldRG.pgm";
    //TestReadImage(input5, output);
    //output = "output//noisy_region.pgmRG.pgm";
    //TestReadImage(input4, output);


    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage, * avgimage;
    image = ReadPNMImage(filename);

    //Task1
    //outimage = OTSU(image);
	//avgimage = AverageFilter(image, 5);
	//outimage = OTSU(avgimage);

    //Task2
	//outimage = PartialOtusImage(image);
	//outimage = Partition(image,image->Height/2,image->Width/3);
    //outimage = LoG(image);


    //Task3
    //outimage = MAverThresholdingImage(image);

    //Task4
    //outimage =RegionGrowing(imagke, 70);

    SavePNMImage(outimage, outfilename);
    return(0);
}



double abs(double a, double b) {
    if (a < b)
        return b - a;
    else
        return a - b;
}

void Findhistogram(int image->Height, int image->Width, Image* image, long histogram[]) {
	short k;
	for (int i = 0; i < image->Height; i++) {
		for (int j = 0; j < image->Width; j++) {
			k = image->data[i * image->Width + j];
			histogram[k] = histogram[k] + 1;
		}
	}
}

void Findpi(int image->Height, int image->Width, long histogram[], double pi[]) {
	for (int i = 0; i < 255; ++i) {
		pi[i] = (double)histogram[i] / (double)(image->Height * image->Width);
	}
}

double p1(int k, double pi[]) {
	double sum = 0.0;
	for (int i = 0; i <= k; i++) {
		sum += pi[i];
	}
	return sum;
}

double m(int k, double pi[]) {
	double sum = 0.0;
	for (int i = 0; i <= k; i++) {
		sum += i * pi[i];
	}
	return sum;
}

double sigma(int k, double mg, double pi[]) {
	double p1k = p1(k, pi);
	double mk = m(k, pi);
	if (p1k < (1e-10) || (1 - p1k) < (1e-10))
		return 0.0;
	else
		return pow(mg * p1k - mk, 2) / (p1k * (1 - p1k));
}

Image* Otsu(Image* image) {
	int size, count = 0;
	long histogram[256];
	double pi[256],var[256];
	double mg;
	short k = 0;
	double max = 0.0;
	double k_star;

    Image* outimage = CreateNewImage(image, (char*)"Otus Image", image->Width, image->Height);

    unsigned char* tempin = image->data;
    unsigned char* tempout = outimage->data;

    Findhistogram(image->Height, image->Width, image, histogram);
    Findpi(image->Height, image->Width, histogram, pi);
	mg = m(255, pi);

	for (int i = 0; i < 256; i++) {
		var[i] = var(i, mg, pi);
	}

	max = var[0];
	count = 1;
	k = 0;
	for (int i = 0; i < 256; i++) {
		if (abs(var[i], max) < 1e-10) {
			k += i;
			count++;
		}
		else if (var[i] > max) {
			max = var[i];
			count = 1;
			k = i;
		}
	}
	k_star = k / count;
	printf("%1f\n", k_star);
	for (int i = 0; i < image->Height; i++) {
		for (int j = 0; j < image->Width; j++) {
			if (tempin[i * image->Width + j] < k_star)
				tempout[i * image->Width + j] = 0;
			else
				tempout[i * image->Width + j] = 255;
		}
	}
	return(outimage);
}

Image* AverageFilter(Image* image, int size) {
	Image* outImage;
	int i, j, p, q, temp;
	unsigned char* tempin, * tempout;

	outImage = CreateNewImage(image, "#AverageFilter", image->Width, image->Height);
	tempout = outImage->data;
	tempin = image->data;

	for (i = size / 2; i < image->Height - size / 2; i++) {
		for (j = size / 2; j < image->Width - size / 2; j++) {
			temp = 0;
			for (p = -1; p < size - 1; p++) {
				for (q = -1; q < size - 1; q++) {
					temp = temp + tempin[(i + p) * image->Width + (j + q)];
				}
			}
			tempout[i * outImage->Width + j] = temp / (size * size);
		}
	}
	return outImage;
}

int Median(int num[], int m) {
	int i, j;
	int temp;
	for (i = 0; i < m / 2 + 1; i++) {
		for (j = i + 1; j < m; j++) {
			if (num[j] < num[i]) {
				temp = num[j];
				num[j] = num[i];
				num[i] = temp;
			}
		}
	}
	return num[m / 2];
}

Image* PartialOtus(Image* image) {
    int size = image->Height * image->Width;
    Image* outimage = CreateNewImage(image, (char*)"#PartialOtus Image", image->Width ,image->Height );
    unsigned char* tempin = image->data;
    unsigned char* tempout = outimage->data;
    int* p, * p1, * ptemp;
    double* m1, * m2,* var;
    int min, max, i, j, m, n, k, mg;
    int areaW = image->Width / 3;
    int areaH = image->Height / 2;

    int areaS = areaW * areaH;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++) {
            * ptemp = (int*)malloc(areaS * sizeof(int));
            for (m = 0; m < areaH; m++) {
                for ( n = 0; n < areaW; n++) {
                    ptemp[m * areaW + n] = tempin[(i * areaH + m) * image->Width + j * areaW + n];
                }
            }
            p = (int*)malloc(sizeof(int) * 256);
            for (i = 0; i < areaS; i++) {
                p[ptemp[i]]++;
            }

            p1 = (int*)malloc(sizeof(int) * 256);
            for (i = 0; i < 256; i++) {
                for (j = 0; j < i; j++) {
                    p1[i] += p[j];
                }
            }

            mg = 0;
            for (i = 0; i < 256; i++) {
                mg += i * p[i];
            }

            m1 = (int*)malloc(sizeof(int) * 256);
            for (int i = 0; i < 256; i++) {
                for (int j = 0; j <= i; j++) {
                    m1[i] += j * p[j];
                }
                if (p1[i] != 0) {
                    m1[i] /= p1[i];
                }
                else {
                    m1[i] = 0;
                }
            }

            m2 = (double*)malloc(sizeof(double) * 256);

            for (int i = 0; i < 256; i++) {
                for (int j = i + 1; j < 256; j++) {
                    m2[i] += j * p[j];
                }
                if ((areaS - p1[i]) != 0) {
                    m2[i] /= (areaS - p1[i]);
                }
                else {
                    m2[i] = 0;
                }
            }

            var = (double*)malloc(sizeof(double) * 256);
            for (i = 0; i < 256; i++) {

                var[i] = p1[i] * pow(m1[i] - mg, 2) + (1 - p1[i]) * pow(m2[i] - mg, 2);
            }

            max = 0, min = 255;
            for (i = 0; i < areaS; i++) {
                if (tempin[i] > max) {
                    max = tempin[i];
                }
                if (tempin[i] < min) {
                    min = tempin[i];
                }
            }

            int k = min + 1;
            for (int i = min + 1; i < max - 1; i++) {
                if (var[i] > var[k]) {
                    k = i;
                }
            }
            for (int m = 0; m < areaH; m++) {
                for (int n = 0; n < areaW; n++) {
                    int temp = tempin[(i * areaH + m) * image->Width + j * areaW + n];
                    if (temp > k) {
                        tempout[(i * areaH + m) * image->Width + j * areaW + n] = 255;
                    }
                    else {
                        tempout[(i * areaH + m) * image->Width + j * areaW + n] = 0;
                    }
                }
            }
            free(var);
            free(p);
            free(p1);
            free(m1);
            free(m2);
            free(ptemp);

        }
    }
    return outimage;
}

Image* MovingAverageThresholding(Image* image) {
    int size = image->Height * image->Width;

    Image* outimage = CreateNewImage(image, (char*)"MovingAverageThresholding", image->Width, image->Height);

    unsigned char* tempin = image->data;
    unsigned char* tempout = outimage->data;

    double diff;
    int n = 20;
    float c = 0.5, m0 = (float)tempin[0] / n, m1;

    for (int i = 0; i < image->Height; i++) {
        for (int j = 0; j < image->Width; j++) {
            diff= 0.0;

            int i * image->Width + j = i * image->Width + j;
            if (i * image->Width + j < n + 1) {
                diff = tempin[i * image->Width + j];
            }
            else {
                diff = tempin[i * image->Width + j] - tempin[i * image->Width + j - n - 1];
            }

            diff *= 1 / n;
            m1 = m0 + diff;
            m0 = m1;

            if (tempin[i * image->Width + j] > round(m1 * c)) {
                tempout[i * image->Width + j] = 255;
            }
            else {
                tempout[i * image->Width + j] = 0;
            }
        }
    }

    return outimage;
}

Image* RegionGrowing(Image* image, float threshold) {
    int size = image->Heigh * image-Width;
    unsigned char* tempin = image->data;
    unsigned char* tempout = outimage->data;
    Image* outimage = CreateNewImage(image, image->Heigh, image-Width, (char*)"#RGrowing Image");

    typedef struct PTS {
        int pix;
        int x;
        int y;
        int pnum;
    } PTS;

    stack<PTS> seed;
    vector<PTS> poly;
    PTS top;
    PTS* pts = new PTS[image-Width * image->Heigh];
    vector<vector<PTS>> polys;



    for (int i = 0; i < image->Heigh; ++i) {
        for (int j = 0; j < image-Width; ++j) {
            pts[i * image-Width + j].pix = input[i * image-Width + j];
            pts[i * image-Width + j].x = i;
            pts[i * image-Width + j].y = j;
            pts[i * image-Width + j].pnum = -1;
        }
    }

    for (int i = 0; i < image->Heigh; ++i) {
        for (int j = 0; j < image-Width; ++j) {
            if (pts[i * image-Width + j].pix <= 0) continue;
            if (pts[i * image-Width + j].pnum > -1) continue;
            if (seed.empty() == true) seed.push(pts[i * image-Width + j]);

            while (!seed.empty()) {
                top = seed.top();
                seed.pop();
                poly.push_back(top);
                pts[(int)(top.x * image-Width + top.y)].pnum = polys.size();

                for (int X = -1; X <= 1; ++X) {

                    for (int Y = -1; Y <= 1; ++Y) {

                        if (top.x + X < 0 || top.x + X >= image->Heigh || top.y + Y < 0 || top.y + Y >= image-Width) continue;
                        if (pts[((int)top.x + X) * image-Width + ((int)top.y + Y)].pnum > -1) continue;
                        if (abs(pts[((int)top.x + X) * image-Width + ((int)top.y + Y)].pix - top.pix) <= threshold) {
                            seed.push(pts[((int)top.x + X) * image-Width + ((int)top.y + Y)]);
                            pts[((int)top.x + X) * image-Width + ((int)top.y + Y)].pnum = polys.size();
                        }
                    }
                }
            }
            polys.push_back(poly);
            poly.clear();
        }
    }

    for (int i = 0; i < polys.size(); ++i) {
        for (int j = 0; j < polys[i].size(); ++j) {
            tempout[(int)(polys[i][j].x * image-Width + polys[i][j].y)] = polys[i][j].pix;
        }
    }

    return outimage;
}
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
#define PI 3.141592

void SavePNMImage(Image*, char*);
Image* SwapImage(Image*);
Image* ReadPNMImage(char*);
Image* Translation(Image*);
Image* Rotation(Image*);
Image* BilinearInterpolation(Image*, int);
Image* VerticalShear(Image*, int);
Image* HorizontalShear(Image*, int);
Image* AverageFilter(Image*, int);
Image* MedianFilter(Image*, int);
Image* Laplacian(Image*);
bool isEmpty(int);
int TestReadImage(char*, char*);
int Median(int num[], int m);
int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//camera.pgm";
    char* output;

    //output = "output//lenaTrans--.pgm";
    //TestReadImage(input1, output);
    //output = "output//lenaTrans++.pgm";
    //TestReadImage(input1, output);
    //output = "output//lenaTrans+-.pgm";
    //TestReadImage(input1, output);
    //output = "output//lenaTrans-+.pgm";
    //TestReadImage(input1, output);

    //output = "output//cameraTrans--.pgm";
    //TestReadImage(input2, output);
    //output = "output//cameraTrans++.pgm";
    //TestReadImage(input2, output);
    //output = "output//cameraTrans+-.pgm";
    //TestReadImage(input2, output);
    //output = "output//cameraTrans-+.pgm";
    //TestReadImage(input2, output);

    //output = "output//cameraRot45.pgm";

 /*   output = "output//cameraShearH30.pgm";*/
    //output = "output//cameraShearV30.pgm";
    //output = "output//LenaShearH30.pgm";
    //output = "output//LenaShearV30.pgm";
    output = "output//LenaLap.pgm";
    TestReadImage(input1, output);
    //output = "output//cameraAvg5.pgm";
    //TestReadImage(input2, output);

    //output = "output//cameraMed5.pgm";
    //TestReadImage(input2, output);

    //output = "output//LenaAvg5.pgm";
    //TestReadImage(input2, output);

    //output = "output//cameraMed5.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaMed5.pgm";
    //TestReadImage(input1, output);

    //output = "output//cameraMed3.pgm";
    //TestReadImage(input2, output);

    //output = "output//LenaMed3.pgm";
    //TestReadImage(input1, output);

    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage, * outimageP;

    image = ReadPNMImage(filename);
    //outimage = Translation(image, -100, -100);
    //SavePNMImage(outimage, outfilename);
    //outimage = Translation(image, +100, +100);
    //SavePNMImage(outimage, outfilename);
    //outimage = Translation(image, +100, -100);
    //SavePNMImage(outimage, outfilename);
    //outimage = Translation(image, -100, +100);
    //SavePNMImage(outimage, outfilename);
    //outimage = Rotation(image, 45);
    //SavePNMImage(outimage, outfilename);

    //outimage = HorizontalShear(image, 30);
    //outimage = VerticalShear(image, 30);
    outimage = Laplacian(image);
    SavePNMImage(outimage, outfilename);

    //outimage = AverageFilter(image, 5);
    //SavePNMImage(outimage, outfilename);

    //outimage = MedianFilter(image, 5);
    //SavePNMImage(outimage, outfilename);

    //outimage = MedianFilter(image, 3);
    //SavePNMImage(outimage, outfilename);
    return(0);
}

Image* BilinearInterpolation(Image* image,Image* outImage) {
    int i, j, x, y;
    float scrx, scry, u, v;
    float scaleX, scaleY;
    unsigned char* tempin, * tempout;


    tempout = outImage->data;
    tempin = image->data;

    scaleX = (float)(image->Width - 1) / (outImage->Width - 1);
    scaleY = (float)(image->Height - 1) / (outImage->Height - 1);
    printf("%f", scaleX);

    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            scrx = i * scaleY;
            scry = j * scaleX;
            x = floor(scrx);
            y = floor(scry);
            u = scrx - x;
            v = scry - y;
            tempout[i * outImage->Width + j] = (float)(1 - u) * (float)(1 - v) * (float)tempin[x * image->Width + y]
                + (float)u * (float)(1 - v) * (float)tempin[(x + 1) * image->Width + y]
                + (float)(1 - u) * (float)v * (float)tempin[x * image->Width + y + 1]
                + (float)u * (float)v * (float)tempin[(x + 1) * image->Width + y + 1];
        }
    }
    return outImage;
}

Image* Translation(Image* image, int tx, int ty ) {
    int i, j, sx= tx, sy=ty;
    int size;
    unsigned char* tempin, * tempout;

    Image* outImage;
    outImage = CreateNewImage(image, "#Translation",image->Width+abs(tx), image->Height+abs(ty));
    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            tempout[i * outImage->Width + j] = 0;
        }
    }

    if (ty <= 0) sy = 0;
    if (tx <= 0) sx = 0;

    for (i = sy; i < sy + image->Height; i++) {
        for (j = sx; j < sx + image->Width; j++) {
            tempout[i * outImage->Width + j] = tempin[(i-sy)* image->Width + (j-sx)];
           
        }
    }
    return outImage;
}

Image* Rotation(Image* image, int degree) {
    int i, j, x = 0, y = 0, p, q,temp=0,sum;
    int size= image->Height* image->Width;
    float srx, sry, rx, ry;
    int tempArr[26];
    unsigned char* tempin, * tempout;
    double angle = degree * PI / 180;
    double co = cos(angle);
    double si = sin(angle);

    int rotateX, rotateY;

    rotateX = ceil(image->Width * fabs(si) + image->Height * fabs(co));
    rotateY = ceil(image->Width * fabs(co) + image->Height * fabs(si));

    Image* outImage;
    outImage = CreateNewImage(image, "#NearesNeighbour", rotateX, rotateY);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            tempout[i * outImage->Width + j] = 0;
            
        }
    }

    for (i = 0; i <image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            x = round(i * si  - j * co+(float)(image->Width*co));
            y = round(i * co + j * si);
            tempout[ x * outImage->Width + y] = tempin[i * image->Width + j];
        }
    }
    for (i = 1; i < outImage->Height-1; i++) {
        for (j = 1; j < outImage->Width-1; j++) {
            if (tempout[i * outImage->Width + j] == 0) {
                sum = 0;
                temp = 0;
                for (p = -1; p < 2; p++) {
                    for (q = -1; q < 2; q++) {
                        tempArr[temp++] = tempout[(i + p) * outImage->Width + (j + q)];
                        sum = sum +tempout[(i + p) * outImage->Width + (j + q)];
                        
                    }
                }
                tempout[i * outImage->Width + j] = sum / 8;
            }
        }
    }
    return outImage;
}

Image* HorizontalShear(Image* image, int degree) {
    int i, j, x = 0, y = 0,p,q,temp,sum;
    int tempArr[26];
    unsigned char* tempin, * tempout;
    double angle = degree * PI / 180;
    double ta = tan(angle);
    double co = cos(angle);
    double si = sin(angle);
    int ShearX, ShearY;

    ShearX = ceil(image->Width + image->Height * fabs(ta));
    ShearY = image->Height;

    Image* outImage;
    outImage = CreateNewImage(image, "#NearesNeighbour", ShearX, ShearY);
    tempout = outImage->data;
    tempin = image->data;


    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            tempout[i * outImage->Width + j] = 0;
        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            x = i;
            y = round((-ta) * i + j+(image->Width*ta));
            tempout[x * outImage->Width + y] = tempin[i * image->Width + j];
        }
    }

    for (i = 1; i < outImage->Height - 1; i++) {
        for (j = 1; j < outImage->Width - 1; j++) {
            if (tempout[i * outImage->Width + j] == 0) {
                sum = 0;
                temp = 0;
                for (p = -1; p < 2; p++) {
                    for (q = -1; q < 2; q++) {
                        tempArr[temp++] = tempout[(i + p) * outImage->Width + (j + q)];
                        sum = sum + tempout[(i + p) * outImage->Width + (j + q)];

                    }
                }
                tempout[i * outImage->Width + j] = sum / 8;
            }
        }
    }
    return outImage;
}

Image* VerticalShear(Image* image, int degree) {
    int i, j, x = 0, y = 0,sum,temp,q,p;
    int tempArr[26];
    unsigned char* tempin, * tempout;
    double angle = degree * PI / 180;
    double ta = tan(angle);
    double co = cos(angle);
    double si = sin(angle);
    int ShearX, ShearY;

    ShearX = image->Width;
    ShearY = ceil(image->Height+ image->Width * fabs(ta));;

    Image* outImage;
    outImage = CreateNewImage(image, "#VerticalShear", ShearX, ShearY);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            tempout[i * outImage->Width + j] = 0;
        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            x = round(i-ta*j+image->Width*ta);
            y = j;
            tempout[x * outImage->Width + y] = tempin[i * image->Width + j];
        }
    }
    for (i = 1; i < outImage->Height - 1; i++) {
        for (j = 1; j < outImage->Width - 1; j++) {
            if (tempout[i * outImage->Width + j] == 0) {
                sum = 0;
                temp = 0;
                for (p = -1; p < 2; p++) {
                    for (q = -1; q < 2; q++) {
                        tempArr[temp++] = tempout[(i + p) * outImage->Width + (j + q)];
                        sum = sum + tempout[(i + p) * outImage->Width + (j + q)];

                    }
                }
                tempout[i * outImage->Width + j] = sum / 8;
            }
        }
    }
    return outImage;
}

Image* AverageFilter(Image* image, int size) {
    Image* outImage;
    int i, j, p, q,temp;
    unsigned char* tempin, * tempout;

    outImage = CreateNewImage(image, "#AverageFilter", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;
    
    for (i = size / 2; i < image->Width -size / 2; i++) {
        for (j = size / 2; j < image->Height - size / 2; j++) {
            temp = 0;
            for (p = -1; p < size-1; p++) {
                for (q = -1; q < size-1; q++) {
                    temp = temp+tempin[(i + p) * image->Width + (j + q)];
                }
            }
            tempout[i * outImage->Width + j] = temp /(size * size);
        }
    }
    return outImage;
}

Image* MedianFilter(Image* image, int size) {
    Image* outImage;
    int i, j, p, q,count=0, temp,tempArr[26];
    unsigned char* tempin, * tempout;
    

    outImage = CreateNewImage(image, "#MedianFilter", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = size / 2; i < image->Width - size / 2; i++) {
        for (j = size / 2; j < image->Height - size / 2; j++) {
            temp = 0;
            for (p = -1; p < size - 1; p++) {
                for (q = -1; q < size - 1; q++) {
                    tempArr[temp++]=tempin[(i + p) * image->Width + (j + q)];
                }
            }
            tempout[i * outImage->Width + j] = Median(tempArr,size*size);      
        }
    }
    return outImage;
}

int Median(int num[], int m) {
    int i, j;
    int temp;
    for (i = 0; i < m/2+1; i++) {
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

bool isEmpty(int num[]) {
    int i;
    for (i = 0; i < 9; i++) {
        if (num[i] != 0) return false;
    }
    return true;
}

Image* Laplacian(Image* image) {
    Image* outImage;
    int i, j, p, q, temp;
    unsigned char* tempin, * tempout;

    outImage = CreateNewImage(image, "#Laplacian", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 1; i < image->Width - 1; i++) {
        for (j = 1; j < image->Height - 1; j++) {
            tempout[i * outImage->Width + j] = 4 * tempin[i * outImage->Width + j]
                - tempin[(i + 1) * outImage->Width + j]
                - tempin[(i - 1) * outImage->Width + j]
                - tempin[i * outImage->Width + j + 1]
                - tempin[i * outImage->Width + j - 1];
        }
    }
    return outImage;
}
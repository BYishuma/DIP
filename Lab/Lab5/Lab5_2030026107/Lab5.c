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
Image* DFTSepatate(Image*, int*, int*);
Image* Display(Image*);
Image* IDFTPhaseAngle(Image*, float*, float*);
Image* IDFTMagnitude(Image*, float*, float*);
int TestReadImage(char*, char*);


int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//camera.pgm";
    char* output;

    //output = "output//CameraIFFTMag.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaIFFTMag.pgm";
    //TestReadImage(input1, output);
    //output = "output//CameraDFT.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaDFT.pgm";
    //TestReadImage(input1, output);
    output = "output//CameraIFFTPhase.pgm";
    TestReadImage(input2, output);
    output = "output//LenaIFFTPhase.pgm";
    TestReadImage(input1, output);

    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage,*FFTimage;
    image = ReadPNMImage(filename);

    float* Rg = (float*)calloc(image->Width * image->Height, sizeof(float));
    float* Ig = (float*)calloc(image->Width * image->Height, sizeof(float));
    
    ////FFTimage = DFTSepatate(image);
    outimage = DFTSepatate(image, Rg, Ig);
    outimage= IDFTPhaseAngle(image, Rg, Ig);
    //outimage= IDFTMagnitude(image, Rg, Ig);
    //outimage = Display(FFTimage);
    //outimage = DFTSepatate(image);
    SavePNMImage(outimage, outfilename);

    return(0);
}

Image* DFTSepatate(Image* image, float* Rg, float* Ig) {
    unsigned char* tempin, * tempout;
    int size, i, j, m, n;
    float real, imaginary, sum;
    float theta, res;
    Image* outimage;
    float* Rt = (float*)calloc(image->Width * image->Height, sizeof(float));
    float* It = (float*)calloc(image->Width * image->Height, sizeof(float));
    outimage = CreateNewImage(image, "#Discrete Fourier Transform",image->Width,image->Height);
    tempin = image->data;
    tempout = outimage->data;

    for (i = 0; i < image->Height; i++) {  // Shift image
        for (j = 0; j < image->Height; j++) {
            tempin[image->Width * i + j] = pow((-1), (i + j)) * tempin[image->Width * i + j];
        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            for (m = 0; m < image->Height; m++) {
                theta = (float)(-2) * PI * ((float)i * (float)m) / (float)(image->Height);
                Rt[i* image->Width  + j] += (float)tempin[m*image->Width  + j] * cos(theta);
                It[i* image->Width  + j] += (float)tempin[m*image->Width  + j] * sin(theta);
            }
        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            for (n = 0; n < image->Width; n++) {
                theta = (float)(-2) * PI * ((float)j * (float)n) / (float)image->Width;
                Rg[i *image->Width  + j] += Rt[image->Width * i + n] * cos(theta) - It[image->Width * i + n] * sin(theta);
                Ig[i *image->Width  + j] += Rt[image->Width * i + n] * sin(theta) + It[image->Width * i + n] * cos(theta);
            }
            res = sqrt(pow(Rg[i * image->Width + j], 2) + pow(Ig[i * image->Width + j], 2)) / sqrt(image->Width*image->Height);
            tempout[image->Width * i + j] = res;
        }
    }

    free(Rt);
    free(It);

    return(outimage);
}

Image* IDFTPhaseAngle(Image* image, float* Rg, float* Ig) {
    unsigned char* tempin, * tempout,*temp;
    int size, i, j, m, n;
    float res;
    float real, imaginary, sum;
    float theta,mag;
    Image* outimage, *tempimage;
    outimage = CreateNewImage(image, "#Discrete Fourier Transform", image->Width, image->Height);
    tempout = outimage->data;

    float* Rt = (float*)calloc(image->Width * image->Height, sizeof(float));
    float* It = (float*)calloc(image->Width * image->Height, sizeof(float));

    int imsize = image->Height * image->Width;
    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            mag= sqrt(pow(Ig[i * image->Width + j], 2) + pow(Rg[i * image->Width + j], 2));
            Rg[i * image->Width + j] = Rg[i * image->Width + j] / mag;
            Ig[i * image->Width + j] = Ig[i * image->Width + j] / mag;
        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            for (m = 0; m < image->Height; m++) {
                theta = 2 * PI * (float)((float)i * (float)m / (float)image->Height);
                Rt[i * image->Width + j] += Rg[m * image->Width + j] * cos(theta) - Ig[m * image->Width + j] * sin(theta);
                It[i * image->Width + j] += Rg[m * image->Width + j] * sin(theta) + Ig[m * image->Width + j] * cos(theta);
            }
        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            real = 0;
            imaginary = 0;
            for (n = 0; n < image->Width; n++) {
                theta = 2 * PI * (float)((float)j * (float)n / (float)image->Width);
                real+= Rt[i * image->Width + n] * cos(theta) - It[i * image->Width + n] * sin(theta);
                imaginary += Rt[i * image->Width + n] * sin(theta) + It[i * image->Width + n] * cos(theta);
            }
            res = sqrt(pow(imaginary, 2) + pow(real, 2));
            tempout[image->Width * i + j] = res;

        }
    }

    free(Rt);
    free(It);
    return(outimage);
}

Image* IDFTMagnitude(Image* image, float* Rg, float* Ig) {
    unsigned char* tempin, * tempout, * temp;
    int size, i, j, m, n;
    float res;
    float real, imaginary, sum;
    float theta;
    Image* outimage, * tempimage;
    outimage = CreateNewImage(image, "#Discrete Fourier Transform", image->Width, image->Height);
    tempout = outimage->data;

    float* Rt = (float*)calloc(image->Width * image->Height, sizeof(float));
    float* It = (float*)calloc(image->Width * image->Height, sizeof(float));
    float* mag= (float*)calloc(image->Width * image->Height, sizeof(float));
    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            mag[i * image->Width + j] = sqrt(pow(Ig[i * image->Width + j], 2) + pow(Rg[i * image->Width + j], 2))/(float)(image->Width*image->Height);
            Rg[i * image->Width + j] = mag[i * image->Width + j];
            Ig[i * image->Width + j] = 0;
        
        }
    }
    for (i = 0; i < image->Height; i++) {  // Shift image
        for (j = 0; j < image->Width; j++) {
            Rg[i * image->Width + j] = Rg[i * image->Width + j] * pow((-1), (i + j));
        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            for (m = 0; m < image->Height; m++) {
                theta = 2 * PI * ((float)i * (float)m / (float)image->Height);
                Rt[i * image->Width + j] += Rg[m * image->Width + j] * cos(theta);
                It[i * image->Width + j] += Rg[m * image->Width + j] * sin(theta);
            }

        }
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            real = 0;
            imaginary = 0;
            for (n = 0; n < image->Width; n++) {
                theta = 2 * PI * (float)((float)j * (float)n / (float)image->Width);
                real += Rt[i * image->Width + n] * cos(theta) - It[i * image->Width + n] * sin(theta);
                imaginary += Rt[i * image->Width + n] * sin(theta) + It[i * image->Width + n] * cos(theta);
            }
            res = sqrt(pow(imaginary, 2) + pow(real, 2));
            tempout[image->Width * i + j] = res;
        }
    }
    free(Rt);
    free(It);
    return(outimage);
}


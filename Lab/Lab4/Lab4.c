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
Image* GammaCorrection(Image*, float);
Image* Sobel(Image*);
Image* Laplacian(Image*);
Image* GlobalHistogramEnhancement(Image*);
Image* LocalHistogramEnhancement(Image*);
int TestReadImage(char*, char*);

int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//camera.pgm";
    char* output;

    /*output = "output//LenaSob.pgm"; */
    output = "output//CameraSob.pgm";
    TestReadImage(input2, output);
    output = "output//LenaSob.pgm";
    TestReadImage(input1, output);

    //output = "output//Cameragamma0.1.pgm";
    //TestReadImage(input2, output);
    //output = "output//Lenagamma0.1.pgm";
    //TestReadImage(input1, output);
    //output = "output//Cameragamma0.4.pgm";
    //TestReadImage(input2, output);
    //output = "output//Lenagamma0.4.pgm";
    //TestReadImage(input1, output);
    //output = "output//Cameragamma0.7.pgm";
    //TestReadImage(input2, output);
    //output = "output//Lenagamma0.7.pgm";
    //TestReadImage(input1, output);
    //output = "output//Cameragamma1.pgm";
    //TestReadImage(input2, output);
    //output = "output//Lenagamma1.pgm";
    //TestReadImage(input1, output);
    //output = "output//CameraGlobalhist.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaGlobalhist.pgm";
    //TestReadImage(input1, output);
    //output = "output//CameraLocalhist.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaLocalhist.pgm";
    //TestReadImage(input1, output);
    //output = "output//CameraLap.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaLap.pgm";
    //TestReadImage(input1, output);

    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage, * outimageP;

    image = ReadPNMImage(filename);

    //outimage = Laplacian(image);
    outimage = Sobel(image);
    //outimage = GammaCorrection(image, 0.1);
    //outimage = GammaCorrection(image, 0.4);
    //outimage = GammaCorrection(image, 0.7);
    //outimage = GammaCorrection(image, 1);
    //outimage = GlobalHistogramEnhancement(image);
    //outimage = LocalHistogramEnhancement(image);
    SavePNMImage(outimage, outfilename);

    return(0);
}


Image* Laplacian(Image* image) {
    Image* outImage;
    int i, j, p, q, temp,count;
    unsigned char* tempin, * tempout;

    outImage = CreateNewImage(image, "#Laplacian", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 1; i < image->Width - 1; i++) {
        for (j = 1; j < image->Height - 1; j++) {
            temp = -4 * tempin[i * image->Width + j]
                + tempin[(i - 1) * image->Width + j]
                + tempin[(i + 1) * image->Width + j]
                + tempin[i * image->Width + j - 1]
                + tempin[i * image->Width + j + 1];
            temp = (temp < 0) ? 0: ((temp > 255) ? 255 : temp);
            tempout[i * outImage->Width + j] = temp;
        }
    }
    return outImage;
}

Image* Sobel(Image* image) {
    Image* outImage;
    int i, j, p, q, temp, count, x, y;
    unsigned char* tempin, * tempout;

    outImage = CreateNewImage(image, "#Sobel", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 1; i < image->Width - 1; i++) {
        for (j = 1; j < image->Height - 1; j++) {
            x = (-1) * tempin[(i - 1) * image->Width + (j - 1)]
                + tempin[(i - 1) * image->Width + (j + 1)]
                - 2 * tempin[(i)*image->Width + (j - 1)]
                + 2 * tempin[(i)*image->Width + (j + 1)]
                - 1 * tempin[(i - 1) * image->Width + (j + 1)]
                + tempin[(i + 1) * image->Width + (j + 1)];
            y = (-1) * tempin[(i - 1) * image->Width + (j - 1)]
                - 2 * tempin[(i - 1) * image->Width + j]
                - tempin[(i - 1) * image->Width + (j + 1)]
                + tempin[(i + 1) * image->Width + (j - 1)]
                + 2 * tempin[(i + 1) * image->Width + j]
                + tempin[(i + 1) * image->Width + (j + 1)];
            temp = round(sqrt(x * x + y * y));
            temp = (temp > 255) ? 255 : temp;
            tempout[i * outImage->Width + j] = temp;
        }
    }
    return outImage;
}

Image* GammaCorrection(Image* image, float gamma) {
    Image* outImage;
    int i, j, p, q, temp = 0, sum = 0;
    double var,avg;
    unsigned char* tempin, * tempout;

    outImage = CreateNewImage(image, "#GammaCorrection", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < image->Width; i++) {
        for (j = 0; j < image->Height; j++) {
            tempout[i * outImage->Width + j] = round((float)255 * pow((float)tempin[i * image->Width + j] / 255, gamma));
            sum += tempout[i * outImage->Width + j];
        }
    }
    avg = (float)sum / (image->Width * image->Height);
    for (i = 0; i < image->Width; i++) {
        for (j = 0; j < image->Height; j++) {
            temp=temp+pow((float)(tempout[i * outImage->Width + j] - avg),2);
        }
    }
    var = temp / (outImage->Width * outImage->Height);
    printf("%f ", var);
    return outImage;

}

Image* GlobalHistogramEnhancement(Image* image) {
    Image* outImage;
    int i, j;
    float accProb=0;
    unsigned char* tempin, * tempout;
    int hist[257] = { 0 };
    float prob[257] = { 0 };
    int s[257] = { 0 };
    outImage = CreateNewImage(image, "#HistogramEnhancement", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    //Calculate the accumulating number of each pixel value
    for (i = 0; i < image->Width; i++) {
        for (j = 0; j < image->Height; j++) {
            hist[tempin[i * image->Width + j]] += 1;
        }
    }
    //Transfer to the probability
    for (i = 0; i < 256; i++) {
        prob[i]=(float)hist[i]/(image->Height*image->Width);
    }
    for (i = 0; i < 256; i++) {
        accProb += prob[i];
        s[i] = 255 * accProb;
    }

    for (i = 0; i < image->Width; i++) {
        for (j = 0; j < image->Height; j++) {
            tempout[i * outImage->Width + j] = s[tempin[i * image->Width + j]];
        }
    }
    return outImage;
}

Image* LocalHistogramEnhancement(Image* image) {
    Image* outImage;
    int i, j, p, q,r, temp = 0,sum = 0,count=0;
    float accProb;
    unsigned char* tempin, * tempout;
    int hist[257] = { 0 };
    float prob[257] = { 0 };
    int s[257] = { 0 };
    outImage = CreateNewImage(image, "#HistogramEnhancement", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < image->Height; i =i+ 4) {
        for (j = 0; j < image->Width; j =j + 4) {
            
            accProb = 0;
            memset(hist, 0, sizeof(hist));
            memset(prob, 0, sizeof(prob));
            memset(s, 0, sizeof(s));
            //Calculate the accumulating number of each pixel value
            for (p = 0; p < 4; p++) {
                for (q = 0; q < 4; q++) {
                    hist[tempin[(i+p) * image->Width + j+q]] += 1;
                }
            }
            for (r = 0; r < 256; r++) {
                prob[r] = (float)hist[r] /16;
            }
            for (r = 0; r < 256; r++) {
                accProb += prob[r];
                s[r] = 255 * accProb;
            }
            for (p = 0; p < 4; p++) {
                for (q = 0; q < 4; q++) {
                    tempout[(i+p) * outImage->Width + j+q] = s[tempin[(i + p) * image->Width + j + q]];
                }
            }
        }
    }
    return outImage;
}

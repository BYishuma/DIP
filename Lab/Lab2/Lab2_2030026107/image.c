#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include<string.h>
#include <ctype.h>
#include "proto.h"
#include <math.h>
int imaArr[300][300];

void SavePNMImage(Image*, char*);
Image* SwapImage(Image*);
Image* ReadPNMImage(char*);
int TestReadImage(char*, char*);
Image* AlternativeLineReduction(Image* image);
Image* PixelReplication(Image* image, int outputSize);
Image* NearestNeighbour(Image* image, int multiple);
Image* BilinearInterpolation(Image* image, int multiple);
Image* NegativeImage(Image* image);
Image* AnySizeAdjustment(Image* image, float ratio);
int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//camera.pgm";
    char* output;

    //output = "output//lenaPix.pgm";
    //TestReadImage(input1, output);

    //output = "output//camaraPix.pgm";
    //TestReadImage(input2, output);

    //output = "output//lenaNearest.pgm";
    //TestReadImage(input1, output);

    //output = "output//camaraNearest.pgm";
    //TestReadImage(input2, output);

    //output = "output//lenaBilinear.pgm";
    //TestReadImage(input1, output);

    //output = "output//camaraBilinear.pgm";
    //TestReadImage(input2, output);

    //output = "output//lenaGray.pgm";
    //TestReadImage(input1, output);

    //output = "output//camaraGray.pgm";
    //TestReadImage(input2, output);

    //output = "output//lenaAnyLarge.pgm";
    //TestReadImage(input1, output);

    //output = "output//camaraAnyLarge.pgm";
    //TestReadImage(input2, output);

    output = "output//lenaAnySmall.pgm";
    TestReadImage(input1, output);

    output = "output//camaraAnySmall.pgm";
    TestReadImage(input2, output);

    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage,*outimageP;
    
    image = ReadPNMImage(filename);
    //outimage = PixelReplication(image, 3);
    //
    //outimage = NearestNeighbour(image, 3);
    //
    //outimage = BilinearInterpolation(image, 3);
    // 
    //outimage = NegativeImage(image);

    //outimage = AnySizeAdjustment(image, 1.2);

    outimage = AnySizeAdjustment(image, 0.8);

    //outimage = AlternativeLineReduction(image);

    SavePNMImage(outimage, outfilename);

    return(0);
}

Image* AlternativeLineReduction(Image* image) {
    Image* outImage;
    int i, j, size;
    unsigned char* tempin, * tempout;

    outImage = (Image*)malloc(sizeof(Image));
    outImage->Width = image->Width * 0.5;
    outImage->Height = image->Height * 0.5;
    outImage->Type = image->Type;
    if (outImage->Type == GRAY) size = outImage->Width * 0.5 * outImage->Height * 0.5;
    if (outImage->Type == COLOR) size = outImage->Width * 0.5 * outImage->Height * 0.5 * 3;
    outImage->data = (unsigned char*)malloc(size);

    tempout = outImage->data;
    tempin = image->data;
    for (i = 0; i < image->Height; i += 2) {
        for (j = 0; j < image->Width; j += 2) {
            tempout[(image->Width / 2) * (i / 2) + j / 2] = tempin[image->Width * i + j];
        }
    }
    return outImage;
}

Image* PixelReplication(Image* image, int multiple) {
    int i, j, p, q;
    int size;
    unsigned char* tempin, * tempout;

    Image* outImage;
    outImage = CreateNewImage(image, "#NearesNeighbour", multiple);

    tempout = outImage->data;
    tempin = image->data;

    int count = 0;
    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            for (p = 0; p < multiple; p++) {
                for (q = 0; q < multiple; q++) {
                    tempout[i * multiple * outImage->Width + j * multiple + p * outImage->Width + q] = tempin[i * image->Width + j];                }
            }
        }
    }
    return outImage;
}

Image* AnySizeAdjustment(Image* image, float ratio) {
    int i, j;
    int size;
    unsigned char* tempin, * tempout;

    Image* outImage;
    outImage = (Image*)malloc(sizeof(Image));
    outImage->Width = image->Width * ratio;
    outImage->Height = image->Height * ratio;
    outImage->Type = image->Type;
    if (outImage->Type == GRAY) size = outImage->Width * outImage->Height;
    if (outImage->Type == COLOR) size = outImage->Width * outImage->Height * 3;
    outImage->data = (unsigned char*)malloc(size);

    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            tempout[i * outImage->Width + j] = tempin[(int)(i/ratio)*image->Width+ (int)(j/ratio)];
        }
    }
    return outImage;
}

Image* NearestNeighbour(Image* image,int multiple) {
    int i, j;
    float x, y;
    unsigned char* tempin, * tempout;

    Image* outImage;
    outImage = CreateNewImage(image, "#NearesNeighbour", multiple);
    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            x = round((float)(i) / multiple);
            y = round((float)(j ) / multiple);
            tempout[outImage->Width * i + j] = tempin[(int)(image->Width * (int)x + (int)y)];
        }   
    }
    return outImage;
}

Image* BilinearInterpolation(Image* image, int multiple) {
    int i, j, x, y;
    float scrx, scry, u, v;
    float scaleX,scaleY;
    unsigned char* tempin, * tempout;

    Image* outImage;
    outImage = CreateNewImage(image, "#BilinearInterpolation", multiple);

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

Image* NegativeImage(Image* image) {
    int i, j;
    unsigned char* tempin, * tempout;

    Image* outImage;
    outImage = CreateNewImage(image, "#NegativeImage",1);

    tempout = outImage->data;
    tempin = image->data;

    for (i = 0; i < outImage->Height; i++) {
        for (j = 0; j < outImage->Width; j++) {
            tempout[outImage->Width * i + j]=255- tempin[image->Width * i + j];
        }
    }
    return outImage;
}
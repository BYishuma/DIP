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


#define M_PI 3.141592

void SavePNMImage(Image*, char*);
Image* SwapImage(Image*);
Image* ReadPNMImage(char*);
Image* Dilation(Image* image, int flag);
Image* Erosion(Image* image, int flag);
Image* Opening(Image* image, int flag);
Image* Closing(Image* image, int flag);
Image* BoundaryExtraction(Image* image);
Image* BoundaryOverlappingOnly(Image* image);
Image* ParticalOverlappingOnly(Image* image);
Image* NoneOverlapping(Image* image);

void CountConnectedPixel(Image* image, char* outputfile);
unsigned char* DilationArr(unsigned char* tempin, int height, int width);
int FindPoint(unsigned char* tempin, Pixel* pixel, int value, int height, int width);
int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//noisy_rectangle.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//connected.pgm";
    char* input3 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//bubbles_on_black_background.pgm";
    char* input4 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//licoln.pgm";
    char* input5 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//noisy_fingerprint.pgm";
    char* input6 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//U.pgm";
    char* input7 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lenaD2.pgm";
    char* input8 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lenaD3.pgm";

    char* output;


    //output = "output//noisy_rectangleD4.pgm";
    //TestReadImage(input1, output);
    //output = "output//noisy_fingerprintD4.pgm";
    //TestReadImage(input5, output);

    //output = "output//noisy_rectangleC3.pgm";
    //TestReadImage(input1, output);
    //output = "output//noisy_fingerprintC3.pgm";
    //TestReadImage(input5, output);

    output = "output//licolnBound.pgm";
    TestReadImage(input4, output);
    output = "output//UBound.pgm";
    TestReadImage(input6, output);

    //Task3
    //output = "output//connectedres.txt";
    //TestReadImage(input2, output);

    //Task4
    //output = "output//bubbles_on_black_backgroundPOO.pgm";
    //TestReadImage(input3, output);
    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage, * noiseImage;
    image = ReadPNMImage(filename);

    //outimage = Dilation(image,2);

    /*outimage = Erosion(image, 2);*/

    //outimage = Opening(image, 2);
    outimage = BoundaryExtraction(image, 1);
    SavePNMImage(outimage, outfilename);

    //Task3
    //CountConnectedPixel(image, outfilename);

    //Task4
    //outimage = BoundaryOverlappingOnly(image);
    //outimage = ParticalOverlappingOnly(image);
    //outimage=  NoneOverlapping(Image* image);


    SavePNMImage(outimage, outfilename);


    return(0);
}

Image* Dilation(Image* image,int flag) {
    //flag=1 3x3 structural element
    //flag=2 4-neighbour structural element

    int i, j, p, q;
    unsigned char* tempin, * tempout;
    Image* outimage, *tempimage;

    
    tempin = image->data;
    outimage = CreateNewImage(image, "Dilation",image->Width , image->Height);
    int *temp = (int*)malloc(sizeof(int) * image->Height * image->Width);
    tempout = outimage->data;
    

    for (i = 0; i < image->Height * image->Width; i++) {
        temp[i] = tempin[i];
    }

    for (i = 1; i < image->Height-1; i++) {
        for (j = 1; j < image->Width-1; j++) {
            if (flag == 1) {
                for (p = -1; p < 2; p++) {
                    for (q = -1; q < 2; q++) {
                        if (tempin[(i + p) * image->Width + (j + q)] == 255) {
                            temp[i * image->Width + j] = 255;
                            break;
                        }
                    }
                }
            }
            else {
                 if (tempin[(i - 1) * image->Width + j ] == 255) {
                    temp[i * image->Width + j] = 255;
                    break;
                 }
                 else if(tempin[(i + 1) * image->Width + j] == 255){
                     temp[i * image->Width + j] = 255;
                     break;
                 }
                 else if (tempin[i * image->Width + (j - 1)] == 255) {
                     temp[i * image->Width + j] = 255;
                     break;
                 }
                 else if (tempin[i * image->Width + (j + 1)] == 255) {
                     temp[i * image->Width + j] = 255;
                     break;
                 }
            }
        }
    }
    for (i = 0; i < image->Height * image->Width; i++) {
        tempout[i] = temp[i];
    }

    return outimage;
}

Image* Erosion(Image* image, int flag) {
    //flag=1 3x3 structural element
    //flag=2 4-neighbour structural element

    int i, j, p, q;
    unsigned char* tempin, * tempout;
    Image* outimage, * tempimage;


    tempin = image->data;
    outimage = CreateNewImage(image, "Erosion", image->Width, image->Height);
    int* temp = (int*)malloc(sizeof(int) * image->Height * image->Width);
    tempout = outimage->data;


    for (i = 0; i < image->Height * image->Width; i++) {
        temp[i] = tempin[i];
    }

    for (i = 1; i < image->Height - 1; i++) {
        for (j = 1; j < image->Width - 1; j++) {
            if (flag == 1) {
                for (p = -1; p < 2; p++) {
                    for (q = -1; q < 2; q++) {
                        if (tempin[(i + p) * image->Width + (j + q)] == 0) {
                            temp[i * image->Width + j] = 0;
                            break;
                        }
                    }
                }
            }
            else {
                if (tempin[(i - 1) * image->Width + j] == 0) {
                    temp[i * image->Width + j] = 0;
                    break;
                }
                else if (tempin[(i + 1) * image->Width + j] == 0) {
                    temp[i * image->Width + j] = 0;
                    break;
                }
                else if (tempin[i * image->Width + (j - 1)] == 0) {
                    temp[i * image->Width + j] = 0;
                    break;
                }
                else if (tempin[i * image->Width + (j + 1)] == 0) {
                    temp[i * image->Width + j] = 0;
                    break;
                }
            }
        }
    }
    for (i = 0; i < image->Height * image->Width; i++) {
        tempout[i] = temp[i];
    }

    return outimage;
}

Image* Opening(Image* image, int flag) {
    Image* outimage;
    outimage = Erosion(image, flag);
    outimage = Dilation(outimage, flag);
    return outimage;
}

Image* Closing(Image* image, int flag) {
    Image* outimage;
    outimage = Dilation(image, flag);
    outimage = Erosion(outimage, flag);
    return outimage;
}

Image* BoundaryExtraction(Image* image) {
    //With 3*3 structural element
    int i;
    unsigned char* tempin, * tempout;
    Image* outimage;
    int* temp = (int*)malloc(sizeof(int) * image->Height * image->Width);

    tempin = image->data;
    outimage = Erosion(image, 1);
    tempout = outimage->data;

    for (i = 0; i < image->Height * image->Width; i++) {
        tempout[i] = tempin[i] - tempout[i];
    }

    return outimage;
}

Image* BoundaryOverlappingOnly(Image* image) {
    int i, j, p, q;
    unsigned char* tempin, * tempout;
    Image* outimage, * tempimage;
    int size = image->Height * image->Width;

    tempin = image->data;
    outimage = CreateNewImage(image, "BoundaryOverlappingOnly", image->Width, image->Height);
    unsigned char* temp = (unsigned char*)malloc(sizeof(unsigned char) * image->Height * image->Width);
    tempout = outimage->data;

    for (i = 0; i < size; i++) {
        temp[i] = 0;
    }

    temp[0] = 255;

    while (true) {
        unsigned char* Xk = (unsigned char*)malloc(sizeof(unsigned char) * image->Height * image->Width);
        memcpy(Xk, temp, size);
        temp = DilationArr(temp, image->Height, image->Width);
        for (i = 0; i < size; i++) {
            if (temp[i] > tempin[i]) {
                temp[i] = tempin[i];
            }
        }
        if (memcmp(Xk, temp, size) == 0){
            break;
        }   
        free(Xk);
    }

    for (i = 0; i < size; i++) {
        tempout[i] = temp[i];
    }
    return outimage;
}

Image* ParticalOverlappingOnly(Image* image) {
    Image* outimage;
    int i, j, count = 0,n=400;
    int size = image->Height * image->Width;
    unsigned char* tempin, * tempout;
    unsigned char* A = (unsigned char*)malloc(sizeof(unsigned char) * size);
    unsigned char* B = (unsigned char*)malloc(sizeof(unsigned char) * size);
    unsigned char* X1 = (unsigned char*)malloc(sizeof(unsigned char) * size);
    unsigned char* X2 = (unsigned char*)malloc(sizeof(unsigned char) * size);
    int* num = (int*)malloc(sizeof(int) * size);

    outimage = BoundaryOverlappingOnly(image);
    tempin= image->data;
    tempout = outimage->data;
    Pixel* pixel = (Pixel*)malloc(sizeof(Pixel) * size);

    for (int i = 0; i < size; i++) {
        tempout[i] = tempin[i] - tempout[i];
    }

    memcpy(A, image->data, size);
    memset(X1, 0, sizeof(unsigned char) * size);
    memset(X2, 0, sizeof(unsigned char) * size);
    memset(B, 0, sizeof(unsigned char) * size);
    memset(num, 0, sizeof(int) * size);

    while (FindPoint(A, pixel, 255, image->Height, image->Width)) {
        X1[pixel[0].y * image->Width + pixel[0].x] = 255;

        while (1) {
            unsigned char* Xk = (unsigned char*)malloc(sizeof(unsigned char) * size);
            memcpy(Xk, X1, size);
            X1 = DilationArr(X1, image->Height, image->Width);
            for (int i = 0; i < size; i++) {
                if (X1[i] > A[i]) {
                    X1[i] = A[i];
                }
            }
            if (memcmp(Xk, X1, size) == 0) {
                break;
            }
            free(Xk);
        }
        num[count] = FindPoint(X1, pixel, 255, image->Height, image->Width);

        if (num[count] > n) {
            for (int i = 0; i < size; i++) {
                if (X2[i] < X1[i]) {
                    X2[i] = X1[i];
                }
            }
        }
        count++;
        for (int i = 0; i < size; i++) {
            if (B[i] < X1[i]) {
                B[i] = X1[i];
            }
        }
        memset(X1, 0, sizeof(unsigned char) * size);
        for (int i = 0; i < size; i++) {
            A[i] -= B[i];
            if (A[i] < 100) {
                A[i] = 0;
            }
        }
    }
    for (int i = 0; i < size; i++) {
        tempout[i] = X2[i];
    }

    free(pixel);
    free(A);
    free(B);
    free(X1);
    free(X2);
    free(num);

    return outimage;
}

Image* NoneOverlapping(Image* image) {
    int size = image->Height * image->Width;
    Image* WithP,* WithB;
    Image* outimage = CreateNewImage(image, "NoneOverlapping", image->Width, image->Height);

    WithB = BoundaryOverlappingOnly(image);
    WithP = ParticalOverlappingOnly(image);
    
    unsigned char* tempin = image->data;
    unsigned char* temp1 = WithB->data;
    unsigned char* temp2 = WithP->data;
    unsigned char* tempout = outimage->data;


    for (int i = 0; i < size; i++) {
        tempout[i] = tempin[i] - temp1[i] - temp2[i];
    }

    return outimage;
}

void CountConnectedPixel(Image* image, char* outputfile) {

    int i, j, count = 0;
    int size = image->Height * image->Width;
    unsigned char* tempin= image->data;

    unsigned char* A = (unsigned char*)malloc(sizeof(unsigned char) * size);
    unsigned char* B = (unsigned char*)malloc(sizeof(unsigned char) * size);
    unsigned char* X1 = (unsigned char*)malloc(sizeof(unsigned char) * size);
    int* num = (int*)malloc(sizeof(int) * size);
    
    Pixel* pixel = (Pixel*)malloc(sizeof(Pixel) * size);
    FILE* fp = fopen(outputfile, "w");

    memcpy(A, image->data, size);
    memset(X1, 0, sizeof(unsigned char) * size);
    memset(B, 0, sizeof(unsigned char) * size);
    memset(num, 0, sizeof(int) * size);

    while (FindPoint(A, pixel, 255 , image->Height, image->Width)) {
        X1[pixel[0].y * image->Width + pixel[0].x] = 255;

        while (1) {
            unsigned char *Xk = (unsigned char*)malloc(sizeof(unsigned char) * size);
            memcpy(Xk, X1, size);
            X1 = DilationArr(X1, image->Height, image->Width);
            for (int i = 0; i < size; i++) {
                if (X1[i] > A[i]) {
                    X1[i] = A[i];
                }
            }
            if (memcmp(Xk, X1, size) == 0) {
                break;
            }
            free(Xk);
        }
        num[count] = FindPoint(X1, pixel , 255, image->Height, image->Width);

        fprintf(fp, "The label of the %dth pixel is %d\n", count, num[count]);


        for (i = 0; i < size; i++) {
            if (B[i] < X1[i]) {
                B[i] = X1[i];
            }
        }
        memset(X1, 0, sizeof(unsigned char) * size);
        for (i = 0; i < size; i++) {
            A[i] -= B[i];
            if (A[i] < 100) {
                A[i] = 0;
            }
        }
        count++;
    }

    fclose(fp);
    free(pixel);
    free(A);
    free(B);
    free(X1);
    free(num);
}

int FindPoint(unsigned char* tempin, Pixel* pixel ,int value, int height, int width) {
    int len = 0;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (tempin[i * width + j] == value) {
                pixel[len].x = j; 
                pixel[len].y = i;
                
                pixel[len].value = tempin[i * width + j];
                len++;
            }
        }
    }

    return len;
}

unsigned char* DilationArr(unsigned char* tempin, int height, int width) {
    int i, j, p, q;
    unsigned char* tempout = (unsigned char*)malloc(sizeof(unsigned char) * height * width);


    for (i = 0; i < height * width; i++) {
        tempout[i] = tempin[i];
    }

    for (i = 1; i < height - 1; i++) {
        for (j = 1; j < width - 1; j++) {
            for (p = -1; p < 2; p++) {
                for (q = -1; q < 2; q++) {
                    if (tempin[(i + p) * width + (j + q)] == 255) {
                        tempout[i * width + j] = 255;
                        break;
                    }
                }
            }
        }
    }
    return tempout;
}

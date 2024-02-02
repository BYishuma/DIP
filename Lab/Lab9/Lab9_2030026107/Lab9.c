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
Image* Roberts(Image* image);
Image* Prewitt(Image* image);
Image* Sobel(Image* image);
Image* Canny(Image* image, int Th, int Tl);
Image* LoG(Image* image);
Image* GlobalThresholding(Image* image, float T0);
void* Threshold(unsigned char* tempin, float threshold, int size);



int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//building_original.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//headCT-Vandy.pgm";
    char* input3 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//noisy_fingerprint.pgm";
    char* input4 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//polymersomes.pgm ";

    char* output;

    //Task1
    //output = "output//headCT-VandyRob.pgm";
    //TestReadImage(input2, output);
    //output = "output//building_originalRob.pgm";
    //TestReadImage(input1, output);
    //output = "output//noisy_fingerprintRob.pgm";
    //TestReadImage(input3, output);

    //output = "output//headCT-VandyPre.pgm";
    //TestReadImage(input2, output);
    //output = "output//building_originalPre.pgm";
    //TestReadImage(input1, output);
    //output = "output//noisy_fingerprintPre.pgm";
    //TestReadImage(input3, output);

    //output = "output//headCT-VandySob.pgm";
    //TestReadImage(input2, output);
    //output = "output//building_originalSob.pgm";
    //TestReadImage(input1, output);
    //output = "output//noisy_fingerprintSob.pgm";
    //TestReadImage(input3, output);

    //Task2
    //output = "output//headCT-VandyCan.pgm";
    //TestReadImage(input2, output);
    //output = "output//noisy_fingerprintCan.pgm";
    //TestReadImage(input3, output);

    //output = "output//headCT-Vandylog.pgm";
    //TestReadImage(input2, output);
    //output = "output//noisy_fingerprintlog.pgm";
    //TestReadImage(input3, output);

    //Task3
    output = "output//polymersomesglo.pgm";
    TestReadImage(input4, output);
    output = "output//noisy_fingerprintglo.pgm";
    TestReadImage(input3, output);

    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage, * noiseImage;
    image = ReadPNMImage(filename);

    //Task1
    //outimage = Roberts(image);
    //outimage = Prewitt(image);
    //outimage = Sobel(image);

    //Task2
    //outimage = Canny(image,120,40);
    //outimage = LoG(image);


    //Task3
    outimage = GlobalThresholding(image, 3);

    SavePNMImage(outimage, outfilename);
    return(0);
}

Image* Roberts(Image* image) {
    Image* outimage;
    int i, j, p, q, sumX, sumY, res;
    unsigned char* tempin, * tempout;
    int size = image->Height * image->Width;

    outimage = CreateNewImage(image, "Roberts", image->Width, image->Height);
    tempout = outimage->data;
    tempin = image->data;

    int X[2][2] = { {-1,0},{0,1} };
    int Y[2][2] = { {0,-1},{1,0} };

    for (i = 1; i < image->Height - 1; i++) {
        for (j = 1; j < image->Width; j++) {
            sumX = 0;
            sumY = 0;
            for (p = 0; p < 2; p++) {
                for (q = 0; q < 2; q++) {
                    sumX += tempin[(i + p) * image->Width + j + q] * X[p][q];
                    sumY += tempin[(i + p) * image->Width + j + q] * X[p][q];
                }
            }
            res = abs(sumX) + abs(sumY);
            if (res > 255) res = 255;
            tempout[i * image->Width + j] = res;
        }
    }
    Threshold(tempout, 0.25, size);

    return outimage;
}

Image* Prewitt(Image* image) {
    Image* outimage;
    int i, j, p, q, sumX, sumY, res;
    unsigned char* tempin, * tempout;
    int size = image->Height * image->Width;

    outimage = CreateNewImage(image, "Prewitt", image->Width, image->Height);
    tempout = outimage->data;
    tempin = image->data;

    int X[3][3] = { {-1, -1, -1}, {0, 0, 0}, {1, 1, 1} };
    int Y[3][3] = { {-1, 0, 1}, {-1, 0, 1}, {-1, 0, 1} };

    for (i = 1; i < image->Height - 1; i++) {
        for (j = 1; j < image->Width; j++) {
            sumX = 0;
            sumY = 0;
            for (p = -1; p < 2; p++) {
                for (q = -1; q < 2; q++) {
                    sumX += tempin[(i + p) * image->Width + j + q] * X[p + 1][q + 1];
                    sumY += tempin[(i + p) * image->Width + j + q] * Y[p + 1][q + 1];
                }
            }
            res = abs(sumX) + abs(sumY);
            if (res > 255) res = 255;
            tempout[i * image->Width + j] = res;
        }
    }
    Threshold(tempout, 0.35, size);

    return outimage;
}

Image* Sobel(Image* image) {
    Image* outimage;
    int i, j, p, q, sumX, sumY, res;
    unsigned char* tempin, * tempout;
    int size = image->Height * image->Width;

    outimage = CreateNewImage(image, "Sobel", image->Width, image->Height);
    tempout = outimage->data;
    tempin = image->data;

    int X[3][3] = { {-1, -2, -1}, {0, 0, 0}, {1, 2, 1} };
    int Y[3][3] = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };

    for (i = 1; i < image->Height - 1; i++) {
        for (j = 1; j < image->Width; j++) {
            sumX = 0;
            sumY = 0;
            for (p = -1; p < 2; p++) {
                for (q = -1; q < 2; q++) {
                    sumX += tempin[(i + p) * image->Width + j + q] * X[p + 1][q + 1];
                    sumY += tempin[(i + p) * image->Width + j + q] * Y[p + 1][q + 1];
                }
            }
            res = abs(sumX) + abs(sumY);
            if (res > 255) res = 255;
            tempout[i * image->Width + j] = res;
        }
    }
    Threshold(tempout, 0.35, size);

    return outimage;

}

Image* Canny(Image* image, int Th, int Tl) {
    Image* outimage;
    int i, j, p, q, res, x, y, sumX, sumY;
    float sigma = 0.5;
    double sum;
    int kernel_size = 1 + 2 * ceil(3 *sigma);
    unsigned char* tempin, * tempout;
    int size = image->Height * image->Width;


    outimage = CreateNewImage(image, "Canny", image->Width, image->Height);
    tempout = outimage->data;
    tempin = image->data;

    double* gausFilter = (double*)malloc(sizeof(double) * kernel_size * kernel_size);
    int kcenter = (kernel_size - 1) / 2;
    int X[3][3] = { {-1, -2, -1}, {0, 0, 0}, {1, 2, 1} };
    int Y[3][3] = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };

    //step1 gaussian filter
    for (i = 0; i < kernel_size; i++) {
        for (j = 0; j < kernel_size; j++) {
            double e = -(double)(pow((i - kcenter - 1), 2) + pow((j - kcenter - 1), 2)) / (2 * sigma * sigma);
            gausFilter[i * kernel_size + j] = (double)1 / (2 * PI * sigma * sigma) * exp(e);
        }
    }

    for (i = kcenter; i < image->Height - kcenter; i++) {
        for (j = kcenter; j < image->Width - kcenter; j++) {
            sum = 0;
            for (p = 0; p < kernel_size; p++) {
                for (q = 0; q < kernel_size; q++) {
                    x = i + p - kcenter;
                    y = j + q - kcenter;
                    sum += tempin[x * image->Width + y] * gausFilter[p * kernel_size + q];
                }
            }
            tempout[i * image->Width + j] = sum;
        }
    }

    //step2 Compute the gradient magnitude and angle images
    double* M = (double*)malloc(sizeof(double) * size);
    double* alpha = (double*)malloc(sizeof(double) * size);
    memset(M, 0, sizeof(double) * size);
    memset(alpha, 0, sizeof(double) * size);

    for (i = 1; i < image->Height - 1; i++) {
        for (j = 1; j < image->Width; j++) {
            sumX = 0;
            sumY = 0;
            for (p = -1; p < 2; p++) {
                for (q = -1; q < 2; q++) {
                    sumX += tempin[(i + p) * image->Width + j + q] * X[p + 1][q + 1];
                    sumY += tempin[(i + p) * image->Width + j + q] * Y[p + 1][q + 1];
                }
            }
            if (sumX > 255) sumX = 255;
            if (sumY > 255) sumY = 255;
            M[i * image->Width + j] = sqrt(pow(sumX, 2) + pow(sumX, 2));
            alpha[i * image->Width + j] = atan2(sumY, sumX) * 180 / PI;


        }
    }

    //Step3 
    double* Gn = (double*)malloc(sizeof(double) * size);
    double* temp = (double*)malloc(sizeof(double) * size);
    double* Gnh = (double*)malloc(sizeof(double) * size);
    double* Gnl = (double*)malloc(sizeof(double) * size);

    memset(Gn, 0, sizeof(double) * size);
    memset(Gnh, 0, sizeof(double) * size);
    memset(Gnl, 0, sizeof(double) * size);
    memset(temp, 0, sizeof(double) * size);


    for (int i = 1; i < image->Height - 1; i++) {
        for (int j = 1; j < image->Width - 1; j++) {
            double angle = alpha[i * image->Width + j];
            double ms = M[i * image->Width + j];
            if (angle > -22.5 && angle < 22.5 || angle < 157.5 && angle > -157.5) {
                if (ms < M[(i - 1) * image->Width + j] || ms < M[(i + 1) * image->Width + j]) {
                    Gn[i * image->Width + j] = 0;
                }
                else {
                    Gn[i * image->Width + j] = ms;
                }
            }
            else if (angle > 67.5 & angle < 112.5 || angle < -67.5 && angle > -112.5) {
                if (ms < M[i * image->Width + j - 1] || ms < M[i * image->Width + j + 1]) {
                    Gn[i * image->Width + j] = 0;
                }
                else {
                    Gn[i * image->Width + j] = ms;
                }
            }
            else if (angle < 157.5 && angle > 112.5 || angle > -22.5 && angle < -67.5) {
                if (ms < M[(i - 1) * image->Width + j + 1] || ms < M[(i + 1) * image->Width + j - 1]) {
                    Gn[i * image->Width + j] = 0;
                }
                else {
                    Gn[i * image->Width + j] = ms;
                }
            }
            else {
                if (ms < M[(i - 1) * image->Width + j - 1] || ms < M[(i + 1) * image->Width + j + 1]) {
                    Gn[i * image->Width + j] = 0;
                }
                else {
                    Gn[i * image->Width + j] = ms;
                }
            }
            if (Gn[i * image->Width + j] > Tl) {
                Gnl[i * image->Width + j] = 1;
                temp[i * image->Width + j] = Gn[i * image->Width + j];
            }
            else {
                Gnl[i * image->Width + j] = 0;
            }
            if (Gn[i * image->Width + j] > Th) {
                Gnh[i * image->Width + j] = 1;
                temp[i * image->Width + j] = Gn[i * image->Width + j];
            }
            else {
                Gnh[i * image->Width + j] = 0;
            }
            Gnl[i * image->Width + j] -= Gnh[i * image->Width + j];
        }
    }

    for (int i = 0; i < image->Height; i++) {
        for (int j = 0; j < image->Width; j++) {
            int p5 = i * image->Width + j;
            //if (Gnh[p5] == 1) {
            //    int p1 = (i - 1) * image->Width + (j - 1);
            //    int p2 = (i - 1) * image->Width + j;
            //    int p3 = (i - 1) * image->Width + (j + 1);
            //    int p4 = i * image->Width + (j - 1);
            //    int p6 = i * image->Width + (j + 1);
            //    int p7 = (i + 1) * image->Width + (j - 1);
            //    int p8 = (i + 1) * image->Width + j;
            //    int p9 = (i + 1) * image->Width + (j + 1);
            //    if (Gnl[p1] == 1) {
            //        Gnl[p1] = 1.5;
            //    }
            //    if (Gnl[p2] == 1) {
            //        Gnl[p2] = 1.5;
            //    }
            //    if (Gnl[p3] == 1) {
            //        Gnl[p3] = 1.5;
            //    }
            //    if (Gnl[p4] == 1) {
            //        Gnl[p4] = 1.5;
            //    }
            //    if (Gnl[p5] == 1) {
            //        Gnl[p4] = 1.5;
            //    }
            //    if (Gnl[p6] == 1) {
            //        Gnl[p6] = 1.5;
            //    }
            //    if (Gnl[p7] == 1) {
            //        Gnl[p7] = 1.5;
            //    }
            //    if (Gnl[p8] == 1) {
            //        Gnl[p8] = 1.5;
            //    }
            //    if (Gnl[p9] == 1) {
            //        Gnl[p9] = 1.5;
            //    }
            //}

            //if (Gnh[p5] == 1) {
            //    for (p = -1; p < 2; p++) {
            //        for (q = -1; q < 2; q++) {
            //            if (Gnl[(i + p) * image->Width + j + q] == 1) {
            //                Gnl[(i + p) * image->Width + j + q] == 3;
            //            }
            //        }
            //    }
            //}
        }
    }

    for (int i = 0; i < size; i++) {
        if (Gnl[i] != 1.5 && Gnh[i] != 1) {
            temp[i] = 0;
        }
        outimage->data[i] = temp[i];
    }

    free(Gn);
    free(Gnh);
    free(Gnl);
    free(temp);
    free(M);
    free(alpha);

    return outimage;
}

Image* LoG(Image* image) {
    Image* outimage;
    int i, j, p, q, res, sumX, sumY;
    float sigma = 0.5;
    double sum;
    int kernel_size = 1 + 2 * ceil(3 * sigma);
    unsigned char* tempin, * tempout;
    int size = image->Height * image->Width;


    outimage = CreateNewImage(image, "LoG", image->Width, image->Height);
    tempout = outimage->data;
    tempin = image->data;

    for (int i = 2; i < image->Height - 2; ++i) {
        for (int j = 2; j < image->Width - 2; ++j) {
            int d = 16 * tempin[i * image->Width + j]
                - tempin[(i - 2) * image->Width + j]
                - tempin[(i - 1) * image->Width + j - 1]
                - 2 * tempin[(i - 1) * image->Width + j]
                - tempin[(i - 1) * image->Width + j + 1]
                - tempin[i * image->Width + j - 2]
                - 2 * tempin[i * image->Width + j - 2]
                - 2 * tempin[i * image->Width + j + 1]
                - tempin[i * image->Width + j + 2]
                - tempin[(i + 1) * image->Width + j - 1]
                - 2 * tempin[(i + 1) * image->Width + j]
                - tempin[(i + 1) * image->Width + j + 1]
                - tempin[(i + 2) * image->Width + j];

            if (d > 255)
                tempout[i * image->Width + j] = 255;
            else
                tempout[i * image->Width + j] = 0;
        }
    }
    return(outimage);
}

Image* GlobalThresholding(Image* image,float T0) {
    int i, j, m, n;
    float newT;
    int size = image->Height * image->Width;
    Image* outimage = CreateNewImage(image, "GlobalThresholding", image->Width, image->Height);


    unsigned char* tempin = image->data;
    unsigned char* tempout = outimage->data;

    int min = 255, max = 0;
    int mu1 = 0, mu2 = 0;
    float T = (min + max) / 2.0;
    for (int i = 0; i < size; i++) {
        if (tempin[i] > max) {
            max = tempin[i];
        }
        if (tempin[i] < min) {
            min = tempin[i];
        }
    }

    while (1) {
        mu1 = 0.0;
        mu2 = 0.0;
        m = 0;
        n = 0;

        for (int i = 0; i < image->Height; i++) {
            for (int j = 0; j < image->Width; j++) {
                if (tempin[i * image->Width + j] <= T) {
                    mu1 += tempin[i * image->Width + j];
                    m++;
                }
                else {
                    mu2 += tempin[i * image->Width + j];
                    n++;
                }
            }
        }
        if (m != 0)
            mu1 /= m;
        if (n != 0)
            mu2 /= n;
        newT = (mu1 + mu2) / 2.0;
        if (abs(T - newT) < T0)  break;
        else T = newT;
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {

            if (tempin[i * image->Width + j] > T) {
                tempout[i * image->Width + j] = 0;
            }
            else {
                tempout[i * image->Width + j] = 255;
            }
        }
    }
    return outimage;
}

void* Threshold(unsigned char* tempin, float threshold, int size) {

    int max = 0;

    for (int i = 0; i < size; i++) {
        if (tempin[i] > max)
            max = tempin[i];
    }

    for (int i = 0; i < size; i++) {
        if (tempin[i] > threshold * max)
            tempin[i] = 255;
        else
            tempin[i] = 0;
    }
}


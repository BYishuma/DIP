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
void FFTShift(struct _complex* src, int height, int width);
Image* HomoFilter(Image* image, float high, float low, float cutoff, float c);
Image* AddSinNoise(Image* image, float c);
Image* IBRF(Image* image, float cutoff, float width);
Image* BBRF(Image* image, float radius, float order, float width);
Image* GBRF(Image* image, float radius, float order, float width);
Image* MedianFilter(Image* image, int size);
Image* AriAverageFilter(Image* image, int size);
Image* GeoAverageFilter(Image* image, int size);
Image* AlphaTriMeanFilter(Image* image, int size);
Image* AdapMedFilter(Image* image, int size, int smax);
int Median(int num[], int m);
int TestReadImage(char*, char*);

void FFT(struct _complex*, struct _complex*, int flag, int height, int width);
void FFTShift(struct _complex*, int height, int width);
void SplitArray(struct _complex* , struct _complex* , int x, int y, int flag, int height, int width);


int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//bridge.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//goldhill.pgm";
    char* input3 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm";
    char* input4 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lenaWithNoise.pgm";
    char* input5 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//cameraWithNoise.pgm";
    char* input6 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lenaD1.pgm";
    char* input7 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lenaD2.pgm";
    char* input8 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lenaD3.pgm";

    char* output;


    //output = "output//CameraILPF5.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaILPF5.pgm";
    //TestReadImage(input1, output);

    //output = "output//CameraBLFP80.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaBLFP80.pgm";
    //TestReadImage(input1, output);

    //output = "output//CameraGLFP80.pgm";
    //TestReadImage(input2, output);
    //output = "output//LenaGLFP80.pgm";
    //TestReadImage(input1, output);


    //output = "output//bridgeHomo.pgm";
    //TestReadImage(input1, output);
    //output = "output//goldhillHomo.pgm";
    //TestReadImage(input2, output);
    output = "output//LenaNoise.pgm";
    //output = "output//LenaIBRF5.pgm";
    //TestReadImage(input3, output);

    //output = "output//LenawithN_Rect.pgm";
    //TestReadImage(input4, output);

    //output = "output//cameraNoiMed5.pgm";

    //output = "output//LenaD1Adp.pgm";
    //output = "output//LenaIBRF.pgm";
    //output = "output//LenaBGBRF.pgm";
    TestReadImage(input8, output);

    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage,*noiseImage;
    image = ReadPNMImage(filename);

    //outimage = HomoFilter(image, 2, 0.8, 180, 2);
    noiseImage = AddSinNoise(image, 20);
    /*outimage = BBRF(image, 60, 1, 5);*/
    outimage = GBRF(image, 60, 5, 30);
    //outimage = IBRF(noiseImage, 60,30);


    //SavePNMImage(noiseImage, outfilename);
    outimage = AriAverageFilter(image, 5);
    SavePNMImage(outimage, "output//LenaD3Ari.pgm");
    outimage = MedianFilter(image, 5);
    SavePNMImage(outimage, "output//LenaD3Med.pgm");
    outimage = GeoAverageFilter(image, 3);
    SavePNMImage(outimage, "output//LenaD3Geo.pgm");
    outimage = AlphaTriMeanFilter(image, 3, 1);
    SavePNMImage(outimage, "output//LenaD3Alp.pgm");
    outimage = AdapMedFilter(image, 3, 9);
    SavePNMImage(outimage, "output//LenaD3Ada.pgm");

    //SavePNMImage(noiseImage, outfilename);

    return(0);
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

int cmp(const void* a, const void* b) {
    return *(int*)a - *(int*)b;
}

Image* HomoFilter(Image* image, float high, float low, float cutoff, float c) {
    int i, j;
    Image* inimage, * outimage;
    float dist;
    float filter;
    struct _complex* in = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);
    struct _complex* out = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);


    for (i = 0; i < image->Height * image->Width; i++) {
        in[i].x = 1.0 * image->data[i];
        in[i].y = 0.0;
    }

    FFT(in, out, 1, image->Height, image->Width);

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            dist = sqrt(pow((i - (double)image->Height / 2), 2) + pow((j - (double)image->Width / 2), 2));
            filter = (high - low) * (1 - exp((-c) * (pow(dist / cutoff, 2)))) + low;
            out[i * image->Width + j].x *= filter;
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#HomoFilter", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }

    return (outimage);
}

Image* AddSinNoise(Image* image,float c) {
    int i, j;
    Image * outimage;
    unsigned char *tempin, * tempout;
    float dist;
    float filter;
    int temp, min=255, max=0;
    float res;


    outimage = CreateNewImage(image, "#AddSinNoise", image->Height, image->Width);
    tempin = image->data;
    tempout = outimage->data;

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            res = tempin[i * image->Width + j] + c * sin(c * i) + c * sin(c * j);
            res = (res < 0) ? 0 : ((res > 255) ? 255 : res);
            tempout[i * image->Width + j] = res;
        }
    }

    return (outimage);
}

Image* IBRF(Image* image, float cutoff, float width) {
    int i, j;
    unsigned char* tempin, * tempout;
    Image* inimage, * outimage;

    tempin = image->data;
    float dist;
    struct _complex* in = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);
    struct _complex* out = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);


    for (i = 0; i < image->Height * image->Width; i++) {
        in[i].x = 1.0 * image->data[i];
        in[i].y = 0.0;
    }

    FFT(in, out, 1, image->Height, image->Width);

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            dist = sqrt(pow((i - (double)image->Height / 2), 2) + pow((j - (double)image->Width / 2), 2));
            if (dist >(cutoff-width/2)&& dist<(cutoff + width / 2)) {
                out[i * image->Width + j].x = 0;
                out[i * image->Width + j].y = 0;
            }
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#IBRF", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }

    return (outimage);
}

Image* BBRF(Image* image, float radius, float order, float width) {
    int i, j;
    unsigned char* tempin, * tempout;
    Image* outimage;
    float dist;
    struct _complex* in = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);
    struct _complex* out = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);

    tempin = image->data;
    for (i = 0; i < image->Height * image->Width; i++) {
        in[i].x = 1.0 * image->data[i];
        in[i].y = 0.0;
    }

    FFT(in, out, 1, image->Height, image->Width);

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            dist = sqrt(pow((i - (double)image->Height / 2), 2) + pow((j - (double)image->Width / 2), 2));
            out[i * image->Width + j].x *= 1 / (1 + pow(dist * width / (pow(dist, 2) - pow(radius, 2)), 2 * order));
            out[i * image->Width + j].y *= 1 / (1 + pow(dist * width / (pow(dist, 2) - pow(radius, 2)), 2 * order));
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#BLPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }
    return (outimage);
}

Image* GBRF(Image* image, float radius, float order, float width) {
    int i, j;
    unsigned char* tempin, * tempout;
    Image* outimage;
    float dist;
    tempin = image->data;

    struct _complex* in = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);
    struct _complex* out = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);

    for (i = 0; i < image->Height * image->Width; i++) {
        in[i].x = 1.0 * image->data[i];
        in[i].y = 0.0;
    }

    FFT(in, out, 1, image->Height, image->Width);

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            dist = sqrt(pow((i - (double)image->Height / 2), 2) + pow((j - (double)image->Width / 2), 2));
            out[i * image->Width + j].x *= 1 - exp(-pow((pow(dist, 2) - pow(radius, 2)) / (dist * width), 2));
            out[i * image->Width + j].y *= 1 - exp(-pow((pow(dist, 2) - pow(radius, 2)) / (dist * width), 2));
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#GLPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }
    return (outimage);
}

Image* MedianFilter(Image* image, int size) {
    Image* outImage;
    int i, j, p, q, count = 0, temp, tempArr[26];
    unsigned char* tempin, * tempout;


    outImage = CreateNewImage(image, "#MedianFilter", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = size / 2; i < image->Width - size / 2; i++) {
        for (j = size / 2; j < image->Height - size / 2; j++) {
            temp = 0;
            for (p = -1; p < size - 1; p++) {
                for (q = -1; q < size - 1; q++) {
                    tempArr[temp++] = tempin[(i + p) * image->Width + (j + q)];
                }
            }
            tempout[i * outImage->Width + j] = Median(tempArr, size * size);
        }
    }
    return outImage;
}

Image* AriAverageFilter(Image* image, int size) {
    Image* outImage;
    int i, j, p, q, temp;
    unsigned char* tempin, * tempout;

    outImage = CreateNewImage(image, "#AriAverageFilter", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = size / 2; i < image->Width - size / 2; i++) {
        for (j = size / 2; j < image->Height - size / 2; j++) {
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

Image* GeoAverageFilter(Image* image, int size) {
    Image* outImage;
    int i, j, p, q;
    unsigned char* tempin, * tempout;
    double temp;
    outImage = CreateNewImage(image, "#GeoAverageFilter", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = size / 2; i < image->Width - size / 2; i++) {
        for (j = size / 2; j < image->Height - size / 2; j++) {
            temp = 1;
            for (p = -1; p < size - 1; p++) {
                for (q = -1; q < size - 1; q++) {
                    temp *= (double)tempin[(i + p) * image->Width + (j + q)];
                }
            }
            tempout[i * outImage->Width + j] = pow(temp ,1.0 / (double)(size * size));
        }
    }

    for (int i = 0; i < image->Height * image->Width; i++) {
        tempout[i] = (tempout[i] < 0) ? 0 : ((tempout[i] > 255) ? 255 : tempout[i]);
    }
    return outImage;
}

Image* AlphaTriMeanFilter(Image* image, int size, float d) {
    Image* outImage;
    int i, j, p, q, temp;
    unsigned char* tempin, * tempout;

    outImage = CreateNewImage(image, "#AlphaTriMeanFilter", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = size / 2; i < image->Width - size / 2; i++) {
        for (j = size / 2; j < image->Height - size / 2; j++) {
            temp = 0;
            for (p = -1; p < size - 1; p++) {
                for (q = -1; q < size - 1; q++) {
                    temp = temp + tempin[(i + p) * image->Width + (j + q)];
                }
            }
            tempout[i * outImage->Width + j] = temp / (size * size-d);
        }
    }
    return outImage;
}

Image* AdapMedFilter(Image* image, int size, int smax) {
    Image* outImage;
    int i, j, p, q, count = 0, temp, tempArr[26], min, max, mid, z, a1, a2, b1, b2, n;
    unsigned char* tempin, * tempout;
    int len = size * size;

    outImage = CreateNewImage(image, "#MedianFilter", image->Width, image->Height);
    tempout = outImage->data;
    tempin = image->data;

    for (i = size / 2; i < image->Width - size / 2; i++) {
        for (j = size / 2; j < image->Height - size / 2; j++) {
            temp = 0;

            for (p = -1; p < size - 1; p++) {
                for (q = -1; q < size - 1; q++) {
                    tempArr[temp++] = tempin[(i + p) * image->Width + (j + q)];
                }
            }
            n = size;
            qsort(tempArr, len, sizeof(int), cmp);
            mid = Median(tempArr, size * size);
            min = tempArr[0];
            max = tempArr[len - 1];
            z = tempin[i * image->Width + j];
            a1 = mid - min;
            a2 = mid - max;
            if (a1 > 0 && a2 < 0) {
                b1 = z - min;
                b2 = z - max;
                if (b1 > 0 && b2 < 0)
                    tempout[i * outImage->Width + j] = z;
                else
                    tempout[i * outImage->Width + j] = mid;
            }
            else {
                n += 2;
                if (n > smax) {
                    tempout[i * image->Width + j] = z;
                    break;
                }
            }
        }
    }
    return outImage;
}

void FFT(struct _complex* in, struct _complex* out, int flag, int height, int width) {

    int y, x, i, u, k, n;

    double theta;
    struct _complex w, a0, a1, t;

    if (flag == -1)
        FFTShift(in, height,width);
    for (y = 0; y < height; y++) {
        SplitArray(&in[y *width + 0], &out[y *width + 0],width, 0, 0, height,width);
        for (i = 0; i < log(1.0 *width) / log(2.0); i++) { 
            n = 2 * pow(2.0, i);                
            for (k = 0; k <width / n; k++) {  
                for (u = 0; u < n / 2; u++) {  
                    theta = -1 * 2 * M_PI * u / n;
                    w.x = cos(theta);
                    w.y = flag * sin(theta); 

                    a0 = out[y *width + k * n + u];
                    a1 = out[y *width + k * n + u + n / 2];

                    t.x = w.x * a1.x - w.y * a1.y;
                    t.y = w.x * a1.y + w.y * a1.x;

                    out[y *width + k * n + u].x = a0.x + t.x;
                    out[y *width + k * n + u].y = a0.y + t.y;
                    out[y *width + k * n + u + n / 2].x = a0.x - t.x;
                    out[y *width + k * n + u + n / 2].y = a0.y - t.y;
                }
            }
        }
        if (flag == 1) 
            for (u = 0; u <width; u++) {
                out[y *width + u].x /=width;
                out[y *width + u].y /=width;
            }
    }

    for (x = 0; x <width; x++) {
        SplitArray(&out[0 *width + x], &out[0 *width + x], 0, height, 1, height,width);
        for (i = 0; i < log(1.0 * height) / log(2.0); i++) { 
            n = 2 * pow(2.0, i);    
            for (k = 0; k < height / n; k++) {   
                for (u = 0; u < n / 2; u++) {  
                    theta = -1* 2 * M_PI * u / n; 
                    w.x = cos(theta);
                    w.y = flag * sin(theta);

                    a0 = out[(k * n + u) *width + x];
                    a1 = out[(k * n + u + n / 2) *width + x];

                    t.x = w.x * a1.x - w.y * a1.y;
                    t.y = w.x * a1.y + w.y * a1.x;

                    out[(k * n + u) *width + x].x = a0.x + t.x;
                    out[(k * n + u) *width + x].y = a0.y + t.y;
                    out[(k * n + u + n / 2) *width + x].x = a0.x - t.x;
                    out[(k * n + u + n / 2) *width + x].y = a0.y - t.y;
                }
            }
        }
        if (flag == -1)
            for (u = 0; u < height; u++) {
                out[u *width + x].x = abs(out[u *width + x].x) / height;
                out[u *width + x].y = abs(out[u *width + x].y) / height;
            }
    }
    if (flag == 1) FFTShift(out, height,width);
}

void FFTShift(struct _complex* in, int height, int width) {
    int i, j, res;
    struct _complex tmp;

    for (i = 0; i < height / 2; i++) {
        for (j = 0; j <width; j++) {
            res = ((i + height / 2) % height) *width + (width / 2 + j) %width;
            tmp = in[i * width + j];
            in[i * width + j] = in[res];
            in[res] = tmp;
        }
    }
}

void SplitArray(struct _complex* in, struct _complex* out, int x, int y, int flag, int height, int width) {
    int i;
    struct _complex* t = (struct _complex*)malloc(sizeof(struct _complex) * (flag == 0 ? x / 2 : y / 2));

    if (flag == 0) {
        if (x <= 1) return;

        for (i = 0; i < x / 2; i++) {
            t[i].x = in[i * 2 + 1].x;
            t[i].y = in[i * 2 + 1].y;

            out[i].x = in[i * 2].x;
            out[i].y = in[i * 2].y;
        }
        for (i = 0; i < x / 2; i++) {
            out[i + x / 2].x = t[i].x;
            out[i + x / 2].y = t[i].y;
        }
        SplitArray(out, out, x / 2, y, flag, height, width);
        SplitArray(out + x / 2, out + x / 2, x / 2, y, flag, height,width);
    }
    else {
        if (y <= 1) return;

        for (i = 0; i < y / 2; i++) {
            t[i].x = in[(i * 2 + 1) *width].x;
            t[i].y = in[(i * 2 + 1) *width].y;

            out[i * width].x = in[(i * 2) * width].x;
            out[i *width].y = in[(i * 2) *width].y;
        }

        for (i = 0; i < y / 2; i++) {
            out[(i + y / 2) *width].x = t[i].x;
            out[(i + y / 2) *width].y = t[i].y;
        }

        SplitArray(out, out, x, y / 2, flag, height,width);
        SplitArray(out +width * y / 2, out +width * y / 2, x, y / 2, flag, height,width);
    }
}

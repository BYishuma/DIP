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
void* DFTSepatate(struct _complex* in, struct _complex* out, int flag, int width, int height);
void FFTShift(struct _complex* src, int height, int width);
Image* ILPF(Image* image, float radius);
Image* IHPF(Image* image, float radius);
Image* BLPF(Image*, float, float);
Image* GLPF(Image*, float, float);
Image* BHPF(Image*, float, float);
Image* GHPF(Image*, float, float);
Image* ThresholdFilter(Image*, int);
int TestReadImage(char*, char*);

void FFT(struct _complex*, struct _complex*, int flag, int height, int width);
void FFTShift(struct _complex*, int height, int width);
void SplitArray(struct _complex* , struct _complex* , int x, int y, int flag, int height, int width);


int main(int argc, char** argv)
{
    char* input1 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm";
    char* input2 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//camera.pgm";
    char* input3 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//fingerprint1.pgm";
    char* input4 = "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//fingerprint2.pgm";
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


    output = "output//fingerprint1IHPF60.pgm";
    TestReadImage(input3, output);
    output = "output//fingerprint2IHPF60.pgm";
    TestReadImage(input4, output);

    //output = "output//fingerprint1BHPF06.pgm";
    //TestReadImage(input3, output);
    //output = "output//fingerprint2BHPF60.pgm";
    //TestReadImage(input4, output);

    //output = "output//fingerprint1GHPF60.pgm";
    //TestReadImage(input3, output);
    //output = "output//fingerprint2GHPF60.pgm";
    //TestReadImage(input4, output);
    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage,*FFTimage;
    image = ReadPNMImage(filename);

    //outimage = ILPF(image, 5);
    //outimage = BLPF(image, 80, 2);
    //outimage = GLPF(image, 80, 2);

    //outimage = IHPF(image, 60);
    //outimage = BHPF(image, 60, 2);
    outimage = GHPF(image, 60, 2);
    outimage = ThresholdFilter(outimage, 20);
    SavePNMImage(outimage, outfilename);

    return(0);
}


Image* ILPF(Image* image, float radius) {
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
            if (dist > radius) {
                out[i * image->Width + j].x = 0;
                out[i * image->Width + j].y = 0;
            }
        
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image,"#ILPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }

    return (outimage);
}

Image* BLPF(Image* image, float radius, float order) {
    int i, j;
    unsigned char* tempin, * tempout;
    Image * outimage;
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
            out[i * image->Width + j].x *=  1 / (1 + pow(dist / radius, 2 * order));
            out[i * image->Width + j].y *=  1 / (1 + pow(dist / radius, 2 * order));
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#BLPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }
    return (outimage);
}

Image* GLPF(Image* image, float radius, float order) {
    int i, j;
    unsigned char* tempin, * tempout;
    Image * outimage;
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
            out[i * image->Width + j].x *= exp((-pow(dist / radius, 2 * order) / (2 * pow(order, 2))));
            out[i * image->Width + j].y *= exp((-pow(dist / radius, 2 * order) / (2 * pow(order, 2))));
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#GLPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }
    return (outimage);
}

Image* IHPF(Image* image, float radius) {
    int i, j;
    unsigned char* tempin, * tempout;
    Image * outimage;
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
            if (dist < radius) {
                out[i * image->Width + j].x = 0;
                out[i * image->Width + j].y = 0;
            }

        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#IHPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }
    return (outimage);
}

Image* BHPF(Image* image, float radius, float order) {
    unsigned char* tempin, * tempout;
    Image* outimage;
    float dist;
    tempin = image->data;

    struct _complex* in = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);
    struct _complex* out = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        in[i].x = 1.0 * image->data[i];
        in[i].y = 0.0;
    }

    FFT(in, out, 1, image->Height, image->Width);

    for (int i = 0; i < image->Height; i++) {
        for (int j = 0; j < image->Width; j++) {
            dist = sqrt(pow((i - (double)image->Height / 2), 2) + pow((j - (double)image->Width / 2), 2));
            out[i * image->Width + j].x *= (1 - 1 / (1 + pow(dist / radius, 2 * order)));
            out[i * image->Width + j].y *= (1 - 1 / (1 + pow(dist / radius, 2 * order)));
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#BLPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }

    return (outimage);
}

Image* GHPF(Image* image, float radius, float order) {
    unsigned char* tempin, * tempout;
    Image* outimage;
    float dist;
    tempin = image->data;

    struct _complex* in = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);
    struct _complex* out = (struct _complex*)malloc(sizeof(struct _complex) * image->Height * image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        in[i].x = 1.0 * image->data[i];
        in[i].y = 0.0;
    }

    FFT(in, out, 1, image->Height, image->Width);

    for (int i = 0; i < image->Height; i++) {
        for (int j = 0; j < image->Width; j++) {
            dist = sqrt(pow((i - (double)image->Height / 2), 2) + pow((j - (double)image->Width / 2), 2));
            out[i * image->Width + j].x *= 1 - exp((-pow(dist / radius, 2 * order) / (2 * pow(order, 2))));
            out[i * image->Width + j].y *= 1 - exp((-pow(dist / radius, 2 * order) / (2 * pow(order, 2))));
        }
    }

    FFT(out, out, -1, image->Height, image->Width);

    outimage = CreateNewImage(image, "#GLPF img", image->Height, image->Width);

    for (int i = 0; i < image->Height * image->Width; i++) {
        outimage->data[i] = (out[i].x < 0) ? 0 : ((out[i].x > 255) ? 255 : out[i].x);
    }
    return (outimage);
}

Image* ThresholdFilter(Image* image, int threshold) {
    unsigned char* tempin, * tempout;
    Image* outimage;
    int i;
    tempin = image->data;
    outimage = CreateNewImage(image, "#GLPF img", image->Height, image->Width);
    tempout = outimage->data;
    for (i = 0; i < image->Width * image->Height; i++) {
        if (tempin[i] > threshold) tempout[i] = 255;
        else tempout[i] = tempin[i];
    }

    return (outimage);

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

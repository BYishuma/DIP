#include <stdio.h>
#include <stdlib.h>
//#include <malloc.h>
#include <ctype.h>
#include <string.h>
#include <iostream>
#include <cmath>
#include "proto.h"
#define PI acos(-1)

//Image* DFT(Image* image);
//Image* DFT(Image* image, int *Rg, int *Ig);
Image* DFT(Image* image, double *Rg, double *Ig, double *phase, double *mag);
Image* IDFT1(Image* image, double *Rg, double *Ig);
Image* ILPF(Image* image, double *Rg, double *Ig, double *ilpfRg, double *ilpfIg, int cut);
Image* BLPF(Image* image, double *Rg, double *Ig, double *blpfRg, double *blpfIg, int cut, int order);
Image* GLPF(Image* image, double *Rg, double *Ig, double *glpfRg, double * glpfIg, int cut);
int TestReadImage(char *filename, char *outfilename);
double Dist(int i, int j, int x, int y);

int main(int argc, char **argv)
{
   
   char* input = "/Users/MSJ/Desktop/Digital Image Processing/PGM_IMAGES/lena.pgm";


   char* output = "...";
   TestReadImage(input, output);
   return(0);
}


int TestReadImage(char *filename, char *outfilename)
{
    
    Image *image;
    Image *outimage;
    double* Rg = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* Ig = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* ilpfRg = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* ilpfIg = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* blpfRg = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* blpfIg = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* glpfRg = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* glpfIg = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* phase = (double*)calloc(image->Width*image->Height, sizeof(double));
    double* mag = (double*)calloc(image->Width*image->Height, sizeof(double));

    image=ReadPNMImage(filename);
   

//------------------------Discrete Fourier Transform------------------------
    outimage = DFT(image, Rg, Ig, phase, mag);
    SavePNMImage(outimage, "/Users/MSJ/Desktop/output/dftLenatest.pgm");
//--------------------------------ILPF--------------------------------------
    outimage = ILPF(image, Rg, Ig, ilpfRg, ilpfIg, 10);  // 10, 60, 150, 500
    SavePNMImage(outimage, "/Users/MSJ/Desktop/output/ilpf10.pgm");
    outimage = IDFT1(image, ilpfRg, ilpfIg);
    SavePNMImage(outimage, "/Users/MSJ/Desktop/output/ilpfLena10.pgm");
//--------------------------------BLPF--------------------------------------
//   outimage =  BLPF(image, Rg, Ig, blpfRg, blpfIg, 250, 2);  // 10, 60, 150, 500
//    SavePNMImage(outimage, "/Users/MSJ/Desktop/output/blp250_2.pgm");
//    outimage = IDFT1(image, blpfRg, blpfIg);
//    SavePNMImage(outimage, "/Users/MSJ/Desktop/output/blpfLena250_2.pgm");
////--------------------------------GLPF--------------------------------------
//    outimage = GLPF(image, Rg, Ig, glpfRg, glpfIg, 250);  // 10, 60, 150, 500
//    SavePNMImage(outimage, "/Users/MSJ/Desktop/output/glpf250.pgm");
//    outimage = IDFT1(image, glpfRg, glpfIg);
//    SavePNMImage(outimage, "/Users/MSJ/Desktop/output/glpfLena250.pgm");



   
    free(Rg);
    free(Ig);
    free(phase);
    free(mag);
    free(ilpfRg);
    free(ilpfIg);
    free(blpfRg);
    free(blpfIg);
    free(glpfRg);
    free(glpfIg);
   return(0);
}



Image *DFT(Image* image, double *Rg, double *Ig, double* phase, double* mag) {
   unsigned char *tempin, *tempout;
   int size, i, j, m, n;
   float theta;
   double res1, real, imaginary, sum;
   Image *outimage;
   outimage = CreateNewImage(image, "#Discrete Fourier Transform");
   tempin = image->data;
   tempout = outimage->data;
   int w = image->Width;
   int h = image->Height;

   double* Rt = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* It = (double*)calloc(image->Width*image->Height, sizeof(double));


   for(i = 0; i < h; i++) {  // Shift image
       for(j = 0; j < h; j++) {
           tempin[w * i + j] = (int)pow((-1),(i + j))*tempin[w * i + j];
       }
   }

   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           for(m = 0; m < h; m++) {
               theta = (-2)*PI*(i*m) / h;
               Rt[w * i + j] += tempin[w * m + j] * cos(theta);
               It[w * i + j] += tempin[w * m + j] * sin(theta);
           }
       }
   }

   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = res1 = 0;
           for(n = 0; n < w; n++) {
               theta = (-2)*PI*(j*n) / w;
               real += Rt[w * i + n] * cos(theta) - It[w * i + n] * sin(theta);
               imaginary += Rt[w * i + n] * sin(theta) + It[w * i + n] * cos(theta);
           }
           Rg[w * i + j] = real;
           Ig[w * i + j] = imaginary;
           
           // Magnitude of DFT
           res1 = sqrt(pow(real,2) + pow(imaginary,2));
           mag[w * i + j] = res1;  // Save the original magnitude
           res1 = res1 / 1000;
           res1 = (res1 < 0) ? 0 : ((res1 > 255) ? 255 : res1);
           tempout[w * i + j] = res1;
           // Phase angle of DFT
           phase[w * i + j] = atan(double(imaginary)/double(real));  // [-pi/2,pi/2]
           if(Rg[w * i + j] < 0 && Ig[w * i + j] > 0 && phase[w * i + j] < 0) {
               phase[w * i + j] += PI;  // Change 4th domain to 2th domain
           } else if (Rg[w * i + j] < 0 && Ig[w * i + j] < 0 && phase[w * i + j] > 0) {
               phase[w * i + j] += PI;  // Change 1th domain to 3th domain
           }
//            printf("%f", cos(phase[w * i + j]));
       }
   }

       free(Rt);
       free(It);

       return(outimage);
   }

double Dist(int i, int j, int x, int y){
    double dist = sqrt(pow((i - x), 2) + pow((j - y), 2));
    return(dist);
}


Image* ILPF(Image* image, double *Rg, double *Ig, double *ilpfRg, double * ilpfIg, int cut) {
    unsigned char *tempin, *tempout;
    int size, i, j, m, n;
    double dist, res;
    Image *outimage;
    outimage = CreateNewImage(image, "#Ideal Lowpass Filter");
    tempin = image->data;
    tempout = outimage->data;
    int w = image->Width;
    int h = image->Height;
    
    // Add filter on Rg and Ig
    for(i = 0; i < h; i++) {
        for(j = 0; j < w; j++) {
            res = 0;
            dist = Dist(i, j, h/2, w/2);
            if(dist <= cut){
                ilpfRg[w * i + j] = 1 * Rg[w * i + j];
                ilpfIg[w * i + j] = 1 * Ig[w * i + j];
            } else {
                ilpfRg[w * i + j] = 0;
                ilpfIg[w * i + j] = 0;
            }
        }
    }
    
    // Sample ILPF
    for(i = 0; i < h; i++) {
        for(j = 0; j < w; j++) {
            res = 0;
            tempin[w * i + j] = 255;
            dist = Dist(i, j, h/2, w/2);
            if(dist <= cut){
                tempout[w * i + j] = 1 * tempin[w * i + j];
            } else {
                tempout[w * i + j] = 0;
            }
        }
    }
            
    return(outimage);
}


Image* BLPF(Image* image, double *Rg, double *Ig, double *blpfRg, double * blpfIg, int cut, int order) {
    unsigned char *tempin, *tempout;
    int size, i, j, m, n;
    double dist, res;
    Image *outimage;
    outimage = CreateNewImage(image, "#Butterworth Lowpass Filter");
    tempin = image->data;
    tempout = outimage->data;
    int w = image->Width;
    int h = image->Height;
    
    // Add filter on Rg and Ig
    for(i = 0; i < h; i++) {
        for(j = 0; j < w; j++) {
            res = 0;
            dist = Dist(i, j, h/2, w/2);
            blpfRg[w * i + j] = (1 / (1+ pow(dist/cut, 2*order))) * Rg[w * i + j];
            blpfIg[w * i + j] = (1 / (1 + pow(dist/cut, 2*order))) * Ig[w * i + j];
        }
    }

    // Sample BLPF
    for(i = 0; i < h; i++) {
        for(j = 0; j < w; j++) {
            res = 0;
            tempin[w * i + j] = 255;  // Make it white
            dist = Dist(i, j, h/2, w/2);
            tempout[w * i + j] = (1 / (1+ pow(dist/cut, 2*order))) * tempin[w * i + j];
        }
    }
    return(outimage);
}


Image* GLPF(Image* image, double *Rg, double *Ig, double *glpfRg, double * glpfIg, int cut) {
    unsigned char *tempin, *tempout;
    int size, i, j, m, n;
    double dist, res;
    Image *outimage;
    outimage = CreateNewImage(image, "#Gaussian Lowpass Filter");
    tempin = image->data;
    tempout = outimage->data;
    int w = image->Width;
    int h = image->Height;
    
    for(i = 0; i < h; i++) {
        for(j = 0; j < w; j++) {
            res = 0;
            dist = Dist(i, j, h/2, w/2);
            
            glpfRg[w * i + j] = exp(-(dist*dist) / (2*cut*cut)) * Rg[w * i + j];
            glpfIg[w * i + j] = exp(-(dist*dist) / (2*cut*cut)) * Ig[w * i + j];

        }
    }
    
    // Sample GLFP
    for(i = 0; i < h; i++) {
        for(j = 0; j < w; j++) {
            res = 0;
            dist = Dist(i, j, h/2, w/2);
            tempin[w * i + j] = 255;  // Make it white
            tempout[w * i + j] = exp(-(dist*dist) / (2*cut*cut)) * tempin[w * i + j];
        }
    }
    return(outimage);
        
}



Image *IDFT1(Image* image, double *Rg, double *Ig) {
   unsigned char *tempin, *tempout;
   int size, i, j, m, n;
   double res1, res2, real, imaginary, sum;
   float theta;
   Image *outimage;
   outimage = CreateNewImage(image, "#Inverse Discrete Fourier Transform");
   tempin = image->data;
   tempout = outimage->data;
   int w = image->Width;
   int h = image->Height;
   
   double* Rt = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* It = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputR = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputI = (double*)calloc(image->Width*image->Height, sizeof(double));
   
   
   // Inverse discrete Fourier Transform
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = 0;
           for(m = 0; m < h; m++) {
               theta = 2*PI*(i*m) / h;
               real += Rg[w * m + j] * cos(theta) - Ig[w * m + j] * sin(theta);
               imaginary += Rg[w * m + j] * sin(theta) + Ig[w * m + j] * cos(theta);
           }
           Rt[w * i + j] = real / h;
           It[w * i + j] = imaginary / h;
//            printf("%d\n", Rt[w * i + j]);
       }
   }

   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = sum = 0;
           for(n = 0; n < w; n++) {
               theta = 2*PI*(j*n) / w;
               real += Rt[w * i + n] * cos(theta) - It[w * i + n] * sin(theta);
               imaginary += Rt[w * i + n] * sin(theta) + It[w * i + n] * cos(theta);
           }
           outputR[w * i + j] = real / w;
           outputI[w * i + j] = imaginary / w;
//            printf("%d\n", outputR[w * i + j]);
           sum = sqrt(pow(real/w,2) + pow(imaginary/w,2));
           sum = (sum < 0) ? 0 : ((sum > 255) ? 255 : sum);
           tempout[w * i + j] = sum ;
       }
   }
    

   for(i = 0; i < h; i++) {  // Shift image
       for(j = 0; j < h; j++) {
           tempout[w * i + j] = tempout[w * i + j] / (int)pow((-1),(i + j));
//            printf("%d\n", tempout[w * i + j]);
       }
   }
    



   free(Rt);
   free(It);
   free(outputR);
   free(outputI);


   return(outimage);
}





Image *IDFT2(Image* image, double *phase, double *mag, double *Rg, double *Ig) {
   unsigned char *tempin, *tempout;
   int size, i, j, m, n;
   double res1, res2, real, imaginary, sum;
   float theta;
   Image *outimage;
   outimage = CreateNewImage(image, "#Inverse Discrete Fourier Transform");
   tempin = image->data;
   tempout = outimage->data;
   int w = image->Width;
   int h = image->Height;
   
   double* Rt = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* It = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputR = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputI = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* R = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* I = (double*)calloc(image->Width*image->Height, sizeof(double));


   
   
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           R[w * i + j] = mag[w * i + j] * cos(phase[w * i + j]);
           I[w * i + j] = mag[w * i + j] * sin(phase[w * i + j]);
           
//            if(R[w * i + j] != Rg[w * i + j]){
//                printf("NO\n");
//            }
           
           printf("%f\n", R[w * i + j]);
           printf("%f\n", Rg[w * i + j]);
//
//            printf("%f\n", I[w * i + j]);
//            printf("%f\n", Ig[w * i + j]);

       }
   }

   

   // Inverse discrete Fourier Transform
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = 0;
           for(m = 0; m < h; m++) {
               theta = 2*PI*(i*m) / h;
               real += R[w * m + j] * cos(theta) - I[w * m + j] * sin(theta);
               imaginary += R[w * m + j] * sin(theta) + I[w * m + j] * cos(theta);
           }
           Rt[w * i + j] = real / h;
           It[w * i + j] = imaginary / h;
//            printf("%d\n", It[w * i + j]);
//            printf("%d\n", Rt[w * i + j]);
       }
   }
   
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = sum = 0;
           for(n = 0; n < w; n++) {
               theta = 2*PI*(j*n) / w;
               real += Rt[w * i + n] * cos(theta) - It[w * i + n] * sin(theta);
               imaginary += Rt[w * i + n] * sin(theta) + It[w * i + n] * cos(theta);
           }
           outputR[w * i + j] = real / w;
           outputI[w * i + j] = imaginary / w;

//            printf("%d\n", outputR[w * i + j]);
           sum = sqrt(pow(real/w,2) + pow(imaginary/w,2));
           sum = (sum < 0) ? 0 : ((sum > 255) ? 255 : sum);
           tempout[w * i + j] = sum;
       }
   }
   
   
   for(i = 0; i < h; i++) {  // Shift image
       for(j = 0; j < h; j++) {
           tempout[w * i + j] = tempout[w * i + j] / (int)pow((-1),(i + j));
//            printf("%d\n", tempout[w * i + j]);
       }
   }
       
   free(Rt);
   free(It);
   free(R);
   free(I);
   free(outputR);
   free(outputI);

   return(outimage);
}



Image *IDFT3(Image* image, double *phase) {
   unsigned char *tempin, *tempout;
   int size, i, j, m, n;
   double res1, res2, real, imaginary, sum;
   float theta;
   Image *outimage;
   outimage = CreateNewImage(image, "#Inverse Discrete Fourier Transform");
   tempin = image->data;
   tempout = outimage->data;
   int w = image->Width;
   int h = image->Height;

   double* blank = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* R = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* I = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* Rt = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* It = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputR = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputI = (double*)calloc(image->Width*image->Height, sizeof(double));
   
   
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           blank[w * i + j] = 255;
       }
   }
   
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           res1 = abs(blank[w * i + j]) * cos(phase[w * i + j]);
           res2 = abs(blank[w * i + j]) * sin(phase[w * i + j]);

           R[w * i + j] = res1;
           I[w * i + j] = res2;

//            printf("%f\n", R[w * i + j]);
//            printf("%f\n", Rg[w * i + j]);
       

       }
   }
   
   
   // Inverse discrete Fourier Transform
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = 0;
           for(m = 0; m < h; m++) {
               theta = 2*PI*(i*m) / h;
//                real += (int)Rg[w * m + j] * cos(theta) - Ig[w * m + j] * sin(theta);
//                imaginary += (int)Rg[w * m + j] * sin(theta) + Ig[w * m + j] * cos(theta);
               real += R[w * m + j] * cos(theta) - I[w * m + j] * sin(theta);
               imaginary += R[w * m + j] * sin(theta) + I[w * m + j] * cos(theta);
           }
           Rt[w * i + j] = real / h;
           It[w * i + j] = imaginary / h;
       }
   }

   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = sum = 0;
           for(n = 0; n < w; n++) {
               theta = 2*PI*(j*n) / w;
               real += Rt[w * i + n] * cos(theta) - It[w * i + n] * sin(theta);
               imaginary += Rt[w * i + n] * sin(theta) + It[w * i + n] * cos(theta);
           }
           outputR[w * i + j] = real / w;
           outputI[w * i + j] = imaginary / w;
           sum = sqrt(pow(real/w,2) + pow(imaginary/w,2));
           sum = (sum < 0) ? 0 : ((sum > 255) ? 255 : sum);
           tempout[w * i + j] = sum;
       }
   }
   
   for(i = 0; i < h; i++) {  // Shift image
       for(j = 0; j < h; j++) {
           tempout[w * i + j] = tempout[w * i + j] / (int)pow((-1),(i + j));
       }
   }
   

   free(Rt);
   free(It);
   free(R);
   free(I);
   free(blank);


   return(outimage);
}


Image *IDFT4(Image* image, double *mag) {
   unsigned char *tempin, *tempout;
   int size, i, j, m, n;
   double res1, res2, real, imaginary, sum;
   float theta;
   Image *outimage;
   outimage = CreateNewImage(image, "#Inverse Discrete Fourier Transform");
   tempin = image->data;
   tempout = outimage->data;
   int w = image->Width;
   int h = image->Height;

   double* blank = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* R = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* I = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* Rt = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* It = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputR = (double*)calloc(image->Width*image->Height, sizeof(double));
   double* outputI = (double*)calloc(image->Width*image->Height, sizeof(double));
   
   
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           blank[w * i + j] = 255;
       }
   }
   
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           res1 = abs(mag[w * i + j]) * cos(blank[w * i + j]);
           res2 = abs(mag[w * i + j]) * sin(blank[w * i + j]);

           R[w * i + j] = res1;
           I[w * i + j] = res2;

//            printf("%f\n", R[w * i + j]);
//            printf("%f\n", Rg[w * i + j]);

       }
   }
   
   
   // Inverse discrete Fourier Transform
   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = 0;
           for(m = 0; m < h; m++) {
               theta = 2*PI*(i*m) / w;
//                real += (int)Rg[w * m + j] * cos(theta) - Ig[w * m + j] * sin(theta);
//                imaginary += (int)Rg[w * m + j] * sin(theta) + Ig[w * m + j] * cos(theta);
               real += R[w * m + j] * cos(theta) - I[w * m + j] * sin(theta);
               imaginary += R[w * m + j] * sin(theta) + I[w * m + j] * cos(theta);
           }
           Rt[w * i + j] = real / h;
           It[w * i + j] = imaginary / h;
       }
   }

   for(i = 0; i < h; i++) {
       for(j = 0; j < w; j++) {
           real = imaginary = sum = 0;
           for(n = 0; n < w; n++) {
               theta = 2*PI*(j*n) / h;
               real += Rt[w * i + n] * cos(theta) - It[w * i + n] * sin(theta);
               imaginary += Rt[w * i + n] * sin(theta) + It[w * i + n] * cos(theta);
           }
           outputR[w * i + j] = real / w;
           outputI[w * i + j] = imaginary / w;
           sum = sqrt(pow(real/w,2) + pow(imaginary/w,2));
           sum = (sum < 0) ? 0 : ((sum > 255) ? 255 : sum);
           tempout[w * i + j] = sum;
       }
   }
   
   for(i = 0; i < h; i++) {  // Shift image
       for(j = 0; j < h; j++) {
           tempout[w * i + j] = tempout[w * i + j] / (int)pow((-1),(i + j));
       }
   }
   

   free(Rt);
   free(It);
   free(R);
   free(I);
   free(blank);


   return(outimage);
}







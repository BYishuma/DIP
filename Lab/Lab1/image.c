#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include<string.h>
#include <ctype.h>
#include "proto.h"
int imaArr[300][300];


void SavePNMImage(Image*, char*);
Image* SwapImage(Image*);
Image* ReadPNMImage(char*);
Image* CreateNewImage(Image*, char* comment);
Image* MedianFilter(Image* image);
Image* AverageFilter(Image* image);
int TestReadImage(char*, char*);
void GetArray(Image* image);
int Median(int num[], int m);
int Average(int num[], int m);
int main(int argc, char** argv)
{
    Image* inimage1, * inimage2, * outimage1, * outimage2, * outimage3, * outimage4;
    //TestReadImage("D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm", "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena2.pgm");
    //SwapImage("D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm", "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena3.pgm");
    inimage1=ReadPNMImage("D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena.pgm");
    inimage2 = ReadPNMImage("D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//noise.pgm");

    //outimage=SwapImage(image);
    outimage1=MedianFilter(inimage1);
    outimage2 = AverageFilter(inimage1);
    outimage3 = MedianFilter(inimage2);
    outimage4 = AverageFilter(inimage2);
    SavePNMImage(outimage1, "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena1.pgm");
    SavePNMImage(outimage2, "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//lena2.pgm");
    SavePNMImage(outimage3, "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//noise1.pgm");
    SavePNMImage(outimage4, "D://大三下//Digital Image Processing//PGM_IMAGES//PGM_IMAGES//noise2.pgm");



    return(0);
}

int TestReadImage(char* filename, char* outfilename)
{
    Image* image;
    Image* outimage;

    image = ReadPNMImage(filename);
    outimage = SwapImage(image);
    SavePNMImage(outimage, outfilename);
    
    return(0);
}

Image* SwapImage(Image* image)
{
    unsigned char* tempin, * tempout;
    int i, size=0;
    Image* outimage;

    outimage = CreateNewImage(image, "#Test SwapImage");
    tempin = image->data;
    tempout = outimage->data;
    if (image->Type == GRAY)   size = image->Width * image->Height;
    else if (image->Type == COLOR) size = image->Width * image->Height * 3;
    

    
    for (i = 0; i < size; i++)
    {
        *tempout = *tempin;
        tempin++;
        tempout++;
    }
    return(outimage);
}

/*******************************************************************************/
//Read PPM image and return an image pointer                                   
/**************************************************************************/
Image * ReadPNMImage(char* filename)
{
    char ch;
    int  maxval, Width, Height;
    int size, num, j;
    FILE* fp;
    Image* image;
    int num_comment_lines = 0;


    image = (Image*)malloc(sizeof(Image));

    if ((fp = fopen(filename, "rb")) == NULL) {
        printf("Cannot open %s\n", filename);
        exit(0);
    }

    printf("Loading %s ...", filename);

    if (fscanf(fp, "P%c\n", &ch) != 1) {
        printf("File is not in ppm/pgm raw format; cannot read\n");
        exit(0);
    }
    if (ch != '6' && ch != '5') {
        printf("File is not in ppm/pgm raw format; cannot read\n");
        exit(0);
    }

    if (ch == '5')image->Type = GRAY;  // Gray (pgm)
    else if (ch == '6')image->Type = COLOR;  //Color (ppm)
    /* skip comments */
    ch = getc(fp);
    j = 0;
    while (ch == '#')
    {
        image->comments[num_comment_lines][j] = ch;
        j++;
        do {
            ch = getc(fp);
            image->comments[num_comment_lines][j] = ch;
            j++;
        } while (ch != '\n');     /* read to the end of the line */
        image->comments[num_comment_lines][j - 1] = '\0';
        j = 0;
        num_comment_lines++;
        ch = getc(fp);            /* thanks, Elliot */
    }

    if (!isdigit((int)ch)) {
        printf("Cannot read header information from ppm file");
        exit(0);
    }

    ungetc(ch, fp);               /* put that digit back */

    /* read the width, height, and maximum value for a pixel */
    fscanf(fp, "%d%d%d\n", &Width, &Height, &maxval);

    /*
    if (maxval != 255){
      printf("image is not true-color (24 bit); read failed");
      exit(0);
    }
    */

    if (image->Type == GRAY)
        size = Width * Height;
    else  if (image->Type == COLOR)
        size = Width * Height * 3;
    image->data = (unsigned char*)malloc(size);
    image->Width = Width;
    image->Height = Height;
    image->num_comment_lines = num_comment_lines;
    //printf(Width);
    //printf(Height);

    if (!image->data) {
        printf("cannot allocate memory for new image");
        exit(0);
    }

    num = fread((void*)image->data, 1, (size_t)size, fp);
    //printf("Complete reading of %d bytes \n", num);
    if (num != size) {
        printf("cannot read image data from file");
        exit(0);
    }

    for(j=0;j<image->num_comment_lines;j++){
          printf("%s\n",image->comments[j]);
          }


    fclose(fp);

    /*-----  Debug  ------*/

    if (image->Type == GRAY)printf("..Image Type PGM\n");
    else printf("..Image Type PPM Color\n");
    
    //printf("Width %d\n", Width);
    //printf("Height %d\n",Height);
    //printf("Size of image %d bytes\n",size);
    //printf("maxvalue %d\n", maxval);
    
    return(image);
}

void SavePNMImage(Image* temp_image, char* filename)
{
    int num, j;
    int size;
    FILE* fp;
    //char comment[100];


    printf("Saving Image %s\n", filename);
    fp = fopen(filename, "wb");
    if (!fp) {
        printf("cannot open file for writing");
        exit(0);
    }

    //strcpy(comment,"#Created by Dr Mohamed N. Ahmed");

    if (temp_image->Type == GRAY) {  // Gray (pgm)
        fprintf(fp, "P5\n");
        size = temp_image->Width * temp_image->Height;
    }
    else  if (temp_image->Type == COLOR) {  // Color (ppm)
        fprintf(fp, "P6\n");
        size = temp_image->Width * temp_image->Height * 3;
    }

    for (j = 0; j < temp_image->num_comment_lines; j++)
        fprintf(fp, "%s\n", temp_image->comments[j]);

    fprintf(fp, "%d %d\n%d\n", temp_image->Width, temp_image->Height, 255);

    num = fwrite((void*)temp_image->data, 1, (size_t)size, fp);

    if (num != size) {
        printf("cannot write image data to file");
        exit(0);
    }

    fclose(fp);
}

/*************************************************************************/
/*Create a New Image with same dimensions as input image                 */
/*************************************************************************/

Image* CreateNewImage(Image* image, char* comment)
{
    Image* outimage;
    int size, j;

    outimage = (Image*)malloc(sizeof(Image));

    outimage->Type = image->Type;
    if (outimage->Type == GRAY)   size = image->Width * image->Height;
    else if (outimage->Type == COLOR) size = image->Width * image->Height * 3;

    outimage->Width = image->Width;
    outimage->Height = image->Height;
    outimage->num_comment_lines = image->num_comment_lines;

    /*--------------------------------------------------------*/
    /* Copy Comments for Original Image      */
    for (j = 0; j < outimage->num_comment_lines; j++)
        strcpy(outimage->comments[j], image->comments[j]);

    /*----------- Add New Comment  ---------------------------*/
    strcpy(outimage->comments[outimage->num_comment_lines], comment);
    outimage->num_comment_lines++;


    outimage->data = (unsigned char*)malloc(size);
    if (!outimage->data) {
        printf("cannot allocate memory for new image");
        exit(0);
    }
    return(outimage);
}

void GetArray(Image* image) {
    int i, j;
    int count = 0;
    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            imaArr[i][j] = image->data[count];
            count++;

        }
    }
}

int Median(int num[], int m) {
    int i, j;
    int temp;
    for (i = 0; i < m; i++) {
        for (j = i + 1; j < m; j++) {
            if (num[j] < num[i]) {
                temp = num[i];
                num[i] = num[j];
                num[j] = temp;
            }
        }
    }
    return num[m / 2];
}

int Average(int num[], int m) {
    int sum = 0;
    int i ;
    for (i = 0; i < m; i++) {
        sum = sum + num[i];
    }
    return sum / m;
}

Image* MedianFilter(Image* image) {
    Image* outImage;
    int i, j, k, median;
    int count = 0;
    outImage = SwapImage(image);
    int tempArr[9];
    int tempData[300][300];

    GetArray(image);
    for (i = 1; i < image->Height-1; i++) {
        for (j = 1; j < image->Width - 1; j++) {
            tempArr[0] = imaArr[i - 1][j - 1];
            tempArr[1] = imaArr[i - 1][j];
            tempArr[2] = imaArr[i - 1][j + 1];
            tempArr[3] = imaArr[i][j - 1];
            tempArr[4] = imaArr[i ][j + 1];
            tempArr[5] = imaArr[i + 1][j - 1];
            tempArr[6] = imaArr[i + 1][j];
            tempArr[7] = imaArr[i + 1][j + 1];
            median = Median(tempArr, 8);
            tempData[i][j] = median;
           }
    }
    for (i = 0; i < image->Height; i++) {
        tempData[i][0] = imaArr[i][0];
        tempData[image->Width][0] = imaArr[image->Width][0];
    }
    for (i = 0; i < image->Width; i++) {
        tempData[0][i] = imaArr[0][i];
        tempData[0][image->Height] = imaArr[0][image->Height];
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            outImage->data[count] = tempData[i][j];
            count++;
        }
    }
    return outImage;
}
    
Image* AverageFilter(Image* image) {
    Image* outImage;
    int i, j, k, average;
    int count = 0;
    outImage = SwapImage(image);
    int tempArr[9];
    int tempData[300][300];

    GetArray(image);
    for (i = 1; i < image->Height - 1; i++) {
        for (j = 1; j < image->Width - 1; j++) {
            tempArr[0] = imaArr[i - 1][j - 1];
            tempArr[1] = imaArr[i - 1][j];
            tempArr[2] = imaArr[i - 1][j + 1];
            tempArr[3] = imaArr[i][j - 1];
            tempArr[4] = imaArr[i][j + 1];
            tempArr[5] = imaArr[i + 1][j - 1];
            tempArr[6] = imaArr[i + 1][j];
            tempArr[7] = imaArr[i + 1][j + 1];
            average = Average(tempArr, 8);
            tempData[i][j] = average;
        }
    }
    for (i = 0; i < image->Height; i++) {
        tempData[i][0] = imaArr[i][0];
        tempData[image->Width][0] = imaArr[image->Width][0];
    }
    for (i = 0; i < image->Width; i++) {
        tempData[0][i] = imaArr[0][i];
        tempData[0][image->Height] = imaArr[0][image->Height];
    }

    for (i = 0; i < image->Height; i++) {
        for (j = 0; j < image->Width; j++) {
            outImage->data[count] = tempData[i][j];
            count++;
        }
    }
    return outImage;
}



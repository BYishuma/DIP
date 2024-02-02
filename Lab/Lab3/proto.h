#pragma once

#define GRAY    1
#define COLOR   2


typedef struct Image
{
    int Width;
    int Height;
    int Type;  // Gray=1, Color=2;
    unsigned char* data;
    char comments[15][100];
    int num_comment_lines;
}Image;

Image* CreateNewImage(Image* image, char* comment, int x, int y)
{
    Image* outimage;
    int size, j;

    outimage = (Image*)malloc(sizeof(Image));

    outimage->Type = image->Type;
    if (outimage->Type == GRAY)   size = x * y;
    else if (outimage->Type == COLOR) size = x * y * 3;

    outimage->data = (unsigned char*)malloc(size);
    outimage->Width = x;
    outimage->Height = y;
    outimage->num_comment_lines = image->num_comment_lines;

    /*--------------------------------------------------------*/
    /* Copy Comments for Original Image      */
    for (j = 0; j < outimage->num_comment_lines; j++)
        strcpy(outimage->comments[j], image->comments[j]);

    /*----------- Add New Comment  ---------------------------*/
    strcpy(outimage->comments[outimage->num_comment_lines], comment);
    outimage->num_comment_lines++;


    if (!outimage->data) {
        printf("cannot allocate memory for new image");
        exit(0);
    }
    return(outimage);
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
//Image* SwapImage(Image* image)
//{
//    unsigned char* tempin, * tempout;
//    int i, size = 0;
//    Image* outimage;
//
//    outimage = CreateNewImage(image, "#Test SwapImage", 1);
//    tempin = image->data;
//    tempout = outimage->data;
//    if (image->Type == GRAY)   size = image->Width * image->Height;
//    else if (image->Type == COLOR) size = image->Width * image->Height * 3;
//
//
//
//    for (i = 0; i < size; i++)
//    {
//        *tempout = *tempin;
//        tempin++;
//        tempout++;
//    }
//    return(outimage);
//}

/*******************************************************************************/
//Read PPM image and return an image pointer                                   
/**************************************************************************/
Image* ReadPNMImage(char* filename)
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

    for (j = 0; j < image->num_comment_lines; j++) {
        printf("%s\n", image->comments[j]);
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




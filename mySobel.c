/*
John Yang
Fall 2018 for Robot Vision CAP 4453

This program is a modified version of the sobel algorithm
to output sobel magnitude image and two additional images showing a low and high
threshold applied to the magnitude.

        *** Notes on running the program **

There should be 6 Command line arguments:
4 filenames for input image, magnitude output,
lo threshold output, and hi threshold output, respectively.
The last two arguments should be 40 and 110 as the lo and hi thresholds
in order to match the samples provided by Lobo.


*/
                                    /* Sobel.c */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

        int pic[256][256];
        int outpicx[256][256];
        int outpicy[256][256];
        int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
        int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
        double ival[256][256], loThresh[256][256], hiThresh[256][256], maxival;

int main(int argc, char **argv)
{
        int i,j,p,q,mr,sum1,sum2;
        double loThreshold, hiThreshold;
        FILE *fo1, *fo2, *fo3, *fp1, *fopen();
        char *foobar;
        char *buffer;
        int bufsize = 30;

        // open the original pgm file to be read from.
        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");

        // open the magnitude file to write to.
	      argc--; argv++;
	      foobar = *argv;
	      fo1=fopen(foobar,"wb");

        // open the low threshold file to write to.
	      argc--; argv++;
	      foobar = *argv;
	      fo2=fopen(foobar,"wb");

        // open the high threshold file to write to.
        argc--; argv++;
	      foobar = *argv;
	      fo3=fopen(foobar,"wb");

        // store the low threshold to be used.
        argc--; argv++;
	      foobar = *argv;
	      loThreshold = atof(foobar);

        // store the high threshold to be used.
        argc--; argv++;
	      foobar = *argv;
	      hiThreshold = atof(foobar);

        // scan through and ignore the first 3 header lines input pgm file.
        for(i = 0; i < 3; i++)
        {
          fgets(buffer, bufsize, fp1);
        }

        // read in pixel info from originial image and store in pic table.
        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                  pic[i][j]  &= 0377;
                }
        }

        // perfrom scanning convolution with x and y masks on the original
        // image storing the results in respective outpic tables.
        mr = 1;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             sum1 = 0;
             sum2 = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                   sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
                }
             }
             outpicx[i][j] = sum1;
             outpicy[i][j] = sum2;
          }
        }

        // perform magnitude calculation using input from x and y outpic tables
        // and store result in ival table.  Record max magnitude to be used
        // in scaling later.
        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                      (outpicy[i][j]*outpicy[i][j])));
             if (ival[i][j] > maxival)
                maxival = ival[i][j];

           }
        }

        // write the .pgm header info onto output files.
        fprintf(fo1, "%s\n%s\n%s\n", "P5", "256 256", "255");
        fprintf(fo2, "%s\n%s\n%s\n", "P5", "256 256", "255");
        fprintf(fo3, "%s\n%s\n%s\n", "P5", "256 256", "255");



        // perform scaling procedure on the magnitude table.
        // then apply thresholds to magnitude table to create new
        // tables that reflect those thresholds.
        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
          {
           // Scaling procedure.
           ival[i][j] = (ival[i][j] / maxival) * 255;

           // apply low threshold to magnitude table.
           if(ival[i][j] >= loThreshold) {
             loThresh[i][j] = 255;
           }
           else {
             loThresh[i][j] = 0;
           }

           // apply i threshold to magnitude table.
           if(ival[i][j] >= hiThreshold) {
             hiThresh[i][j] = 255;
           }
           else {
             hiThresh[i][j] = 0;
           }

           // write data to each file.
           fprintf(fo1,"%c",(char)((int)(ival[i][j])));
           fprintf(fo2,"%c",(char)((int)(loThresh[i][j])));
           fprintf(fo3,"%c",(char)((int)(hiThresh[i][j])));


          }
        }

      fclose(fp1);
      fclose(fo1);
      fclose(fo2);
      fclose(fo3);
    /**/
}

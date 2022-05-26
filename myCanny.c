#include <stdio.h>
#include <stdlib.h>              
#include <math.h>
#define  PICSIZE   256
#define  MAXMASK   100
#define  ON        255
#define  OFF       0
#define  HISTOSIZE 255

         int    pic[PICSIZE][PICSIZE];    // original image
         double maskx[MAXMASK][MAXMASK];  // masks to be used in convolution
         double masky[MAXMASK][MAXMASK];
         double convX[PICSIZE][PICSIZE];  // result of convolution with maskX
         double convY[PICSIZE][PICSIZE];  // result of convolution with maskY
         double mag[PICSIZE][PICSIZE];    // image of gradient magnitude.
         double peaks[PICSIZE][PICSIZE];  // image of peaks
         double final[PICSIZE][PICSIZE];  // image after double thresholding
         int histogram[PICSIZE];       // histogram of magnitude values.
         double loPercent = 0.35;

int main(int argc, char **argv)
{
        int     i,j,p,q,mr,centx,centy,moreToDo,areaOfTops;
        double  maskval,exponent,sumX,sumY,sig,maxival,maxval,ZEROTOL, slope;
        double hiThresh, loThresh, percent, cutOff;
        FILE    *fo1, *fo2, *fo3, *fp1, *fopen();
        char    *foobar;

        // open input file to be read.
        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");

        // open output file for magnitudes.
        argc--; argv++;
        foobar = *argv;
        fo1=fopen(foobar,"wb");

        // open output file for peaks.
        argc--; argv++;
        foobar = *argv;
        fo2=fopen(foobar,"wb");

        // open output file for final.
        argc--; argv++;
        foobar = *argv;
        fo3=fopen(foobar,"wb");

        // store the sigma value that will be used for Gaussian smoothing.
        argc--; argv++;
        foobar = *argv;
        sig = atof(foobar);

        // store the percent value that will be used for thresholding.
        argc--; argv++;
        foobar = *argv;
        percent = atof(foobar);

// Output File Prep ===========================================================

        // write pgm header information to output files.
        fprintf(fo1, "%s\n%s\n%s\n", "P5", "256 256", "255");
        fprintf(fo2, "%s\n%s\n%s\n", "P5", "256 256", "255");
        fprintf(fo3, "%s\n%s\n%s\n", "P5", "256 256", "255");

// ============================================================================


        mr = (int)(sig * 3);    // set mask radius accoring to sigma.
        centx = (MAXMASK / 2);  // set centering values for x and y
        centy = (MAXMASK / 2);

// *** find out why the fgets was segfualting for me here but not in sobel.


        // scan through and ignore the first 3 header lines of input pgm file.
        for(i = 0; i < 3; i++)
        {
          while(fgetc(fp1) != '\n')
            ;
        }

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                  pic[i][j]  &= 0377;
                }
        }


// ** double check and make sure i am using p and q in the right places.

        // Compute gaussian values from equation and create the masks.
        for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
              // store the exponent component of the gaussian equation.
              exponent = (exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));

              // compute value for maskx and store in maskx.
              maskval = (q)*exponent;
              (maskx[p+centy][q+centx]) = maskval;

              // compute value for masky and store in masky.
              maskval = (p)*exponent;
              (masky[p+centy][q+centx]) = maskval;
           }
        }

        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sumX = 0;
             sumY = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sumX += pic[i+p][j+q] * maskx[p+centy][q+centx];
                   sumY += pic[i+p][j+q] * masky[p+centy][q+centx];
                }
             }
               convX[i][j] = sumX;
             convY[i][j] = sumY;
          }
        }

        // perform magnitude calculation using input from x and y outpic tables
        // and store result in mag table.  Record max magnitude to be used
        // in scaling later.
        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             mag[i][j]=sqrt((double)((convX[i][j]*convX[i][j]) +
                                      (convY[i][j]*convY[i][j])));
             if (mag[i][j] > maxival)
                maxival = mag[i][j];

           }
        }

        // Scale mag values so they are all within 0 - 255.
        for (i=0;i<256;i++) {
          for (j=0;j<256;j++) {
           mag[i][j] = (mag[i][j] / maxival) * 255;
         }
        }

//================>>>>> Part 2 - Peak Detection <<<<< ==========================

        for(i = mr ; i < 256-mr; i++){
           for(j=mr;j<256-mr;j++){

              // if a value in the convX table is zero, hack it so we can
              // divide by this value safely.
              if((convX[i][j]) == 0.0) {
                 convX[i][j] = .00001;
               }

              // get slope info and do direction checks. We are checking to see
              // in which direction the peak is oriented so that we can run
              // the max test perpendicular to the peak.
              slope = convY[i][j]/convX[i][j];

              if( (slope <= .4142)&&(slope > -.4142)){
                 if((mag[i][j] > mag[i][j-1])&&(mag[i][j] > mag[i][j+1])){
                    peaks[i][j] = 255;
                  }
               }
              else if( (slope <= 2.4142)&&(slope > .4142)) {
                 if((mag[i][j] > mag[i-1][j-1])&&(mag[i][j] > mag[i+1][j+1])) {
                     peaks[i][j] = 255;
                 }
               }
              else if( (slope <= -.4142)&&(slope > -2.4142)) {
                 if((mag[i][j] > mag[i+1][j-1])&&(mag[i][j] > mag[i-1][j+1])) {
                     peaks[i][j] = 255;
                 }
               }
               else {
                 if((mag[i][j] > mag[i-1][j])&&(mag[i][j] > mag[i+1][j])){
                     peaks[i][j] = 255;
                 }
               }
             }
          }

// write peak info to file to show where peak "candidates exist".

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
          {
           fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
          }
        }

//=========== >>>>> Part 3 - Finding Hi and Lo threshold <<<<< =================

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
          {
            histogram[(int)mag[i][j]]++;
          }
        }


        cutOff = percent*(PICSIZE*PICSIZE);
        printf("cutoff: %f\n", cutOff);

        for(i = HISTOSIZE; i >= 1; i--)
        {
          areaOfTops += histogram[i];
          if(areaOfTops > cutOff)
          {
            break;
          }
        }

        hiThresh = i;
        loThresh = loPercent*hiThresh;

        printf("\nHigh Threshold: %03f\nLow Threshold: %03f\n", hiThresh, loThresh);

//============= >>>>> Part 4 - Double Thresholding <<<<< ======================

        // This loop selects all of the highest peaks and includes them in final
        // and eliminates lowest peaks from consideration.
        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
          {
            // if we have a peak and it's greater than threshold, turn it off in
            // peak so it's not revisited and include it in final.
            if(peaks[i][j] == ON)
            {
              if(mag[i][j] > hiThresh)
              {
                peaks[i][j] = OFF;
                final[i][j] = ON;
              }

              // if the magnitude of gradient is lower than low, turn it OFF
              // in peaks so it's not revisited and do not include it as an
              // edge in final.
              else if(mag[i][j] < loThresh) {
                peaks[i][j] = OFF;
                final[i][j] = OFF;
              }
            }
          }
        }

        // Set the moreToDo flag and run while loop.
        // This loop checks for any remaining peaks which, at this point,
        // would be lower than the hi threshold but higher than the lo threshold.
        // If we encounter such a peak, we will check each of its neighbors.
        // If it has a neighbor that is higher than hiThresh and has made it to
        // final image, we will also include it in final.
        // The while loop will continue until no more peaks remain to be checked.
        moreToDo = ON;
        while(moreToDo == ON) {
          moreToDo = OFF;

          for (i=0;i<256;i++) {
            for (j=0;j<256;j++) {
              if (peaks[i][j] == ON)
              {
                for (p = -1; p <= 1; p++) {
                  for (q = -1; q <= 1; q++) {
                    if (final[ i + p ][ j + q ] == ON) {
                        peaks[i][j] = OFF;
                        final[i][j] = ON;
                        moreToDo = ON;
                    }
                  }
                }
              }
            }
          }
        }


//============= Produce Output ======================

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
          {

           // write data to each file.
           fprintf(fo1,"%c",(char)((int)(mag[i][j])));
           fprintf(fo3,"%c",(char)((int)(final[i][j])));
         }
       }

      fclose(fp1);
      fclose(fo1);
      fclose(fo2);
      fclose(fo3);
}

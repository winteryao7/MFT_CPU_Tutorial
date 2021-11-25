/***********************************************
Written by Xiaofeng Meng, Feb 28 2012

This program allows data shift when stacking and output 9times the MAD.

Usage: XmengStackShift File(list of stacking files)      Number of stacking files     Number of data points allowed shift    Name of output file
************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sac.h"
#include "sacio.c"
#include "MaxValue.c"

float   MaxValue (float *array, int index, int n_shift, int size);

static char *outFile = NULL;
FILE        *fp, *NineTimes;
#define     max(a,b) ((a) > (b) ? (a):(b))
#define USAGE "%s [-f list_stacking_sac_file] [-s number_stacking_sac_file]  [-b AllowShift] [-o stacking_sac_file]\n"

int     compare (const void *a, const void *b) {
        const float *da = (const float *) a;
        const float *db = (const float *) b;
        return  (*da > *db) - (*da < *db);
}


int main(int argc, char **argv) {
    
    int StackNumber, NineTimesBreakPoint;

    StackNumber=strtol(argv[2],NULL,10);
    int AllowShift=strtol(argv[3],NULL,10);
    SACHEAD      hd_m[StackNumber], hd_o;
    float        *stack, dt, StackValue, *StackSort, *StackSort2, median, MAD;
    int          DataNumber, i, j;
    
    
    float        *wfcc[StackNumber];
    char         *dataFile[StackNumber][100];
    
    fp=fopen(argv[1],"r");
    if (fp == NULL) {
         printf("Can't open file!\n");
         return -1;
    }

    for (i=0; i<StackNumber; i++) {
        fscanf(fp, "%s", &dataFile[i]);
        fscanf(fp, "%*[^\n]");
    }
    if (feof(fp)) {
       fclose(fp);
    }
    outFile=argv[4];
    for (i=0; i<StackNumber; i++) {
        /*
        printf("Can't\n");
        printf("%s\n",dataFile[i]);
        printf("Can't\n");
        */
        wfcc[i]=malloc(DataNumber*sizeof(float));
        wfcc[i]=read_sac(dataFile[i], &hd_m[i]);
    }

    DataNumber=hd_m[0].npts;
    dt=hd_m[0].delta;

    int          n_shift = 1;

    stack=malloc(DataNumber*sizeof(float));
    StackSort=malloc(DataNumber*sizeof(float));
    StackSort2=malloc(DataNumber*sizeof(float));     
 
    for (i=0; i<DataNumber; i++) {
        for (j=0; j<StackNumber; j++) {
              StackValue=MaxValue(wfcc[j], i, AllowShift, DataNumber);
              stack[i] += StackValue;

            /*
            if (i == 0 || i == DataNumber-1) {
               StackValue=wfcc[j][i];
            }  else {
               StackValue=max(wfcc[j][i-1], wfcc[j][i]);
               StackValue=max(StackValue, wfcc[j][i+1]);
            }   
            stack[i] += StackValue;
            */
        }
    }
    
    for (i=0; i<DataNumber; i++) {
        stack[i]=stack[i]/StackNumber;
        StackSort[i]=stack[i];
    }
    
    qsort(StackSort, DataNumber, sizeof(float), compare);
    if ( DataNumber%2 == 0 ) {
         median=(StackSort[DataNumber/2]+StackSort[DataNumber/2-1])/2;
    }

    if ( DataNumber%2 == 1 ) {
         median=StackSort[(DataNumber-1)/2];
    }    
   
    for (i=0; i<DataNumber; i++) {
         StackSort2[i]=fabsf(StackSort[i]-median);
    }
    
    qsort(StackSort2, DataNumber, sizeof(float), compare);   

    if ( DataNumber%2 == 0 ) {
         MAD=(StackSort2[DataNumber/2]+StackSort2[DataNumber/2-1])/2;
    }

    if ( DataNumber%2 == 1 ) {
         MAD=StackSort2[(DataNumber-1)/2];
    }
    
    float threshold=median+9*MAD;
    /*
    for (i=0; i<DataNumber; i++) {
          if (StackSort[i] >= threshold) {
              NineTimesBreakPoint=i;
              break;
          }
    }
    */

    //NineTimes=fopen(argv[4]".9times","a+");
    printf("mad:%.5f:%.5f\n",median,MAD);
    for (i=0; i<DataNumber; i++) {
         if (stack[i] >= threshold) {
             printf("%.2f %.4f %.2f\n",hd_m[0].b+hd_m[0].delta*i,stack[i],(stack[i]-median)/MAD);
         }
    }
    //fclose(NineTimes);
    int nwindow=rint(DataNumber/n_shift)+1;
    hd_o=sachdr(hd_m[0].delta,nwindow,hd_m[0].b);
    write_sac(outFile,hd_o,stack);
   
   for (j=0; j<StackNumber; j++) {
       free(wfcc[j]);
   }

   free(stack);
   free(StackSort);
   free(StackSort2);
}












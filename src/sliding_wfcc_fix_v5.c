/*********************************************************************
*	sliding_wfcc_fix_v5.c:
*	Compute waveform cross-correlation coeficient between a template
*	and a long continous trace.
*	Usage: src/sliding_wfcc_fix_v5 -f template_sac_file -s long_sac_file -b tBefore_template -a tAfter_template -B tBefore_cont -A tAfter_cont -S -sliding_time_in_sec -O -shift_origin_time -F output -o sacfile [-t -taper] [-D debug] [-F output] [-o sacfile] [-h help]
*
*	Author: Zhigang Peng
*
*	Revision History
*		04/22/07	modified from src_ss.c (by Lupei Zhu), and runDT (by Fenlin Niu)
*		01/02/08	modified from sliding_wfcc.c, this version compute cc without time shift
*	        02/07/08	redo the timing to make sure the output is correct
*		02/07/08	add output option -F 0 (output to screen as two column data) 1 (output to SAC_file)
*		04/04/08	make sure that the output time stamp start with a common for all the trace.
*		09/22/11	fixed a few bugs associated with sliding window as reported by Bogdan Enescu at NIED (bogdan3j@gmail.com)
*		09/25/11	clean up the code so that it can be used for parallel computing
*		09/26/11	add a flag to allow self detection so no need to read the file twice
*		05/03/12        modified by Xiaofeng Meng; change the way of shifting origin time back to version 3.
*		06/06/16        modified by Zhigang Peng; fixed a few minor syntax issues.
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sac.h"
#include "sacio.c"
#include "Complex.h"

#define USAGE "%s -f template_sac_file -s long_sac_file -b tBefore_template -a tAfter_template -B tBefore_cont -A tAfter_cont -S -sliding_time_in_sec -O -shift_origin_time -F output -o sacfile [-t -taper] [-D debug] [-F output] [-o sacfile] [-h help]\n"

// #define TAPER	0.05		/* portion of window being tapered */
/* Note: this taper could be set to 0 so that the raw data is being cross-correlated */

static  char    *dataFile1 = NULL;
static  char    *dataFile2 = NULL;
static  char    *outFile = NULL;

void		_taper_(float *, int, float);
float 	_cal_corr_(int, float *, float *);

int main(int argc, char **argv) {

  extern      int getopt();
  extern      char *optarg;
  extern      int optind;

  int 		i, nn, mm, mt8, t8, error, c, n_shift,ii,nwindow, debug_flag, block=512,n=0;
  int		shift, o_shift_flag, nstart, tmp1, tmp2, output_flag, same_file_flag = 1, taper_flag;
  size_t	wndw_size;
  float		tshift, *output, t0;	/*max. time shift in sec.*/
  float		tBefore_temp, tAfter_temp, arr, max, min, modMaster, modOther,mean_time,tBefore_cont,tAfter_cont, t_shift;
  float 	norm, dt, *src, *master, *other, *trace, *crl, o_shift, tmp_mean, tmp_b1,tmp_b2, *all, taper;
  SACHEAD	hd_m, hd, hd_o;
 
  error = 0;
  tBefore_temp = -5;
  tAfter_temp = 10;
  tshift = 1;
  o_shift_flag = 0;
  debug_flag = 0;
  taper_flag = 0;
  /* input parameters */

    while( (c=getopt( argc, argv, "f:s:b:a:B:A:t:S:O:D:F:o:h" )) != (-1) ) {
        switch( c ) {
        case 'f':
            dataFile1 = optarg;
            break;
        case 's':
            dataFile2 = optarg;
            break;
        case 'b':
            tBefore_temp   = atof(optarg);
            break;
        case 'a':
            tAfter_temp   = atof(optarg);
            break;
        case 'B':
            tBefore_cont   = atof(optarg);
            break;
        case 'A':
            tAfter_cont   = atof(optarg);
            break;
        case 't':
            taper   = atof(optarg);
	    taper_flag = 1;
            break;
        case 'O':
            o_shift   = atof(optarg);
	    o_shift_flag = 1;
            break;
        case 'S':
            t_shift =  atof(optarg);
            break;
        case 'D':
            debug_flag =  atoi(optarg);
            break;
        case 'F':
            output_flag =  atoi(optarg);
            break;
        case 'o':
            outFile =  optarg;
            break;
        case 'h':
            fprintf( stderr, USAGE, argv[0] );
            exit(1);
        default:
            fprintf(stderr, USAGE, argv[0] );
            exit(1);
        }
    }

    if(argc == 1 || dataFile1 == NULL || dataFile2 == NULL) {
        fprintf(stderr, USAGE, argv[0] );
        exit(1);
    }

   if (o_shift_flag == 0) {
      o_shift = 0;
   }
  /* input template trace */
   same_file_flag = strcmp(dataFile1,dataFile2);
  if (same_file_flag == 0) {
	if (debug_flag) printf("template file %s and continuous waveform file %s is the same\n",dataFile1,dataFile2); 
  }
  if ( (src=read_sac(dataFile1,&hd_m)) == NULL ) return -1;
  all = src; // add by zpeng, remember the pointer
  nn = hd_m.npts;
  dt = hd_m.delta;
  n_shift = rint(t_shift/hd_m.delta);

  mm = rint((tAfter_temp-tBefore_temp)/dt);
  wndw_size = mm*sizeof(float);
  if ( debug_flag ) {
    printf("%f %d %d\n",t_shift,mm, n_shift);
  }

  if ( (master=(float *)malloc(wndw_size)) == NULL ||
       (other=(float *)malloc(wndw_size)) == NULL ) {
    fprintf(stderr,"fail to allocation memory for src\n");
    return -1;
  }
  mt8 = rint((tBefore_temp-hd_m.b)/dt);
  tmp1 = rint((tBefore_temp-hd_m.o)/dt);
  if ( debug_flag ) {
  	printf("template waveform start time is %d wndw_size is %d\n",mt8,mm);
  }
// note: this is the true begining time for the template cut data
// so we should add this time shift back
   tmp_b1 = hd_m.b+mt8*hd_m.delta;     //changed by X.Meng

  if ( debug_flag ) {
	printf("template b = %f number from b = %d, number from o = %d\n",hd_m.b,mt8,tmp1);
  }
  if (mt8 < 0) {
    fprintf(stderr,"%s time before arr. is not long enough\n",dataFile1);
    return -1;
  }
  memcpy(master, src+mt8, wndw_size);

   if ( debug_flag ) {
    for(i=0;i<mm;i++) printf("** %d %f %f\n",i,hd_m.b+(mt8+i)*hd_m.delta,master[i]);
   }
// remove mean
  tmp_mean = 0.0;
  for(i=0;i<mm;i++)  tmp_mean+=master[i];
  tmp_mean = tmp_mean/mm;
  for(i=0;i<mm;i++)  master[i]=master[i]-tmp_mean;
  if (taper_flag) _taper_(master,mm,taper);
//  taper(master, mm);

   if ( debug_flag ) {
    for(i=0;i<mm;i++) printf("## %d %f %f\n",i,hd_m.b+(mt8+i)*hd_m.delta,master[i]);
   }
  
//  for(i=0;i<nn;i++) src[i] = 0.;
// read the continuous waveforms
  if (same_file_flag != 0 &&  (trace=read_sac(dataFile2,&hd)) == NULL ) {
    fprintf(stderr,"fail to read data from %s\n",dataFile2);
    return -1;
  }

// set the start and end data point 
// note by zpeng: this part still needs further improvement
//  nstart = rint((tBefore_cont+tmp_b1-hd.b)/hd.delta);
//  if (same_file_flag==0) nstart = rint((tBefore_cont-hd_m.b)/hd_m.delta+mt8); // change tmp_b1 to mt8
//  else nstart = rint((tBefore_cont-hd.b)/hd.delta+mt8);
  if (same_file_flag==0) nstart = rint((tBefore_cont+tmp_b1-hd_m.b)/hd_m.delta); // change tmp_b1 to mt8
  else nstart = rint((tBefore_cont+tmp_b1-hd.b)/hd.delta);  //changed by X.Meng
       printf("start point is %f, %f\n", 1+tmp_b1, hd.b);
/* compute the total sliding window needed */
//  nn = rint((tAfter_cont-tBefore_cont)/t_shift);     //B. Enescu
  if (same_file_flag==0) nn = rint((tAfter_cont-tBefore_cont)/hd_m.delta); //old one
  else nn = rint((tAfter_cont-tBefore_cont)/hd.delta);

  nwindow = rint((nn-mm)/n_shift)+1;
  if (nn < 1) {
    fprintf(stderr,"%s not enough data in %s\n",dataFile2);
    return -1;
  }
  if ( debug_flag ) {
     printf("%f %f %f %f %d %d %d %d\n",tBefore_cont,tAfter_cont,hd_m.delta,hd_m.b,mm,nstart,nn,nwindow);
  }

  if ( output_flag ) {
    n = block;
    output = (float *) calloc(n,sizeof(float));
  }


  for (ii = 1;ii<=nwindow;ii++) {
//    t8 = n_shift*(ii-1)+nstart;
    t8 = n_shift*(ii-1)+nstart; // change by zpeng 09/26/2011
  if ( debug_flag ) {
	printf("--n_shift = %d, ii = %d,nstart = %d, t8 = %d\n",n_shift,ii,nstart,t8);
  }
  if (same_file_flag == 0 ){
//    if (ii == 1) memcpy(other, src-mt8+t8,wndw_size);
//    else memcpy(other, src+t8,wndw_size);
    memcpy(other, all+t8,wndw_size);
  } else {
    memcpy(other, trace+t8,wndw_size);
  }
// remove mean
    tmp_mean = 0.0;
    for(i=0;i<mm;i++)  tmp_mean+=other[i];
    tmp_mean = tmp_mean/mm;
    for(i=0;i<mm;i++)  other[i]=other[i]-tmp_mean;
  if (taper_flag) _taper_(other,mm,taper);
//    taper(other, mm);
   if (ii == 1 && debug_flag ) {
	for(i=0;i<mm;i++) printf("** ii = %d,i = %d, master = %.4f other = %.4f\n",ii,i,master[i],other[i]);
   }
    norm = _cal_corr_(mm,master,other);
// cross-correlation values should be betwen -1 and 1
// this is likely the place that takes most of the CPU time and can be improved

    mean_time = tBefore_cont + (ii-1)*n_shift*dt;  //changed by X.Meng
//    mean_time = tBefore_cont + (ii-1)*n_shift*t_shift - o_shift; //changed by zpeng

// starting time for sac output
   if ( ii == 1 ) t0 = mean_time;
    if ( output_flag == 0 ) {
    	printf("%12.2f %7.4f %d\n",mean_time,norm,ii);
    } else if (output_flag == 1) {
       output[ii-1] = norm;
       if (ii==n) {
          n += block;
          output = realloc(output, n*sizeof(float));
       }
   }
    if ( debug_flag ) {
//     for(i=0;i<mm;i++) printf("** %d %d %f %f %f\n",ii,t8+i,hd_m.b+(t8+i)*dt,master[i],other[i]);
    }
//	printf("## %d %12.2f %7.4f\n",ii,mean_time,norm);
    
    fflush(stdout);
  }  

 if ( output_flag == 1 ) {
//   hd_o = sachdr(dt,nn,t0);
//     hd_o = sachdr(t_shift,nn,t0); //changed by B. Enescu
     hd_o = sachdr(t_shift,nwindow,t0); //changed by Z. Peng
  if ( write_sac(outFile,hd_o,output) == -1 ) {
    fprintf(stderr,"fail to allocation memory for outFile %s\n",outFile);
    return -1;
  }
 }
  free(src); 
  if (same_file_flag != 0) free(trace); 
  free(master); 
  free(other); 
  if ( output_flag ) {
     free(output); 
  }
  return 0;

}

void	_taper_(float *aa, int n, float taper)
{
  int i, m;
  float	tt, pi1;
  m = taper*n;
  pi1 = 3.1415926/m;
  for (i=0; i<m; i++) {
    tt = 0.5*(1.-cos(i*pi1));
    aa[i] *= tt;
    aa[n-i-1] *= tt;
  }
}

float _cal_corr_(n, data, src)
    int     n;
    float   *data, *src;
{
    int     i;
    double  sum, sum1, sum2;

    sum = 0.; sum1=0.; sum2 = 0.;
    for(i = 0; i < n; i++) {
        sum  += (data[i] * src[i]);
        sum1 += (data[i] * data[i]);
        sum2 += (src[i] * src[i]);
    }
//    printf("%f %f %f\n",sum,sum1,sum2);
//    fflush(stdout);
    if((sum1/n) <= 1.e-5 || (sum2/n) <=1.e-5) {
	return(0.0);
    } else {
   	return((float)(sum/sqrt(sum1*sum2)));
    }
}

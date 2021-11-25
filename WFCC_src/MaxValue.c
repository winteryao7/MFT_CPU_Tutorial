#include <stdio.h>
#include <stdlib.h>
#include <string.h>


float   MaxValue (float *array, int index, int n_shift, int size)   {

      int    i;
      float  max;
      
      if (index-n_shift >= 0 && index+n_shift < size) {
            max=array[index-n_shift];

            for (i=index-n_shift; i<index+n_shift; i++) {
                   if (max < array[i+1])   max=array[i+1];
            }
      }
      

      if (index-n_shift < 0) {
            max=array[0];

            for (i=0; i<index+n_shift; i++) {
                   if (max < array[i+1])   max=array[i+1];
            }
      }

     
      if (index+n_shift >= size) {
            max=array[index-n_shift];

            for (i=index-n_shift; i<size-1; i++) {
                   if (max < array[i+1])   max=array[i+1];
            }
      }

      return max;

}

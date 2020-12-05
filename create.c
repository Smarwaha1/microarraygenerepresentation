#include <stdio.h>
#include <math.h>

#include "random.h"

main()
{
  int i,j;
  int no_sample, no_train_sample,no_test_sample1,no_final,no_test_sample;
  int temp, arr[100], count;
  FILE *fp, *fq;

  warmup_random(0.345435);
  
  // Leukemia
  no_sample = 72;
  no_train_sample = 38;
  no_test_sample1 = 14;

  no_test_sample = no_sample - no_train_sample;
  no_final = no_test_sample - no_test_sample1;

  fp = fopen("test_file1.inp","w");
  fq = fopen("test_file2.inp","w");

  for (i=0; i<no_test_sample; i++)
    arr[i] = i;
  for (i=0; i<no_final; i++)
    {
      j = rnd(i,no_test_sample-1);
      
      temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
      
      fprintf(fp,"%d\t",arr[i]+no_train_sample);
    }
  fclose(fp);
  
  for (i=0; i<no_test_sample; i++)
    {
      count=0;
      for (j=0; (j<no_final && i != arr[j]); j++)
	count++;
      if (count >= no_final)
	fprintf(fq,"%d\t",i+no_train_sample);
    }
  fclose(fq);
}

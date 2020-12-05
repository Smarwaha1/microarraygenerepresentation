/*This is the file which initializes the population*/
void init(population *pop_ptr);

void init(population *pop_ptr)
{
  int i, j, k, tmp[optima];
  int max_one;
  

  max_one = 0.1*chrom;
  printf("max_one = %d\n",max_one);
  for(i = 0 ; i < popsize ; i++)
    {
      pop_ptr->ind_ptr = &(pop_ptr->ind[i]);

      for(j = 0; j < chrom; j++)
	pop_ptr->ind_ptr->genes[j] = 0;

      for(j = 0; j < max_one; j++)
	{
	  tmp[j] = rnd(1, chrom);
	  for(k = 0; k < j; k++)
	    {
	      if(tmp[j] == tmp[k])
		{
		  tmp[j] = rnd(1, chrom);
		  k = -1;
		}
	    }
	  pop_ptr->ind_ptr->genes[tmp[j]-1] = 1;
	}
    }
  return;
}

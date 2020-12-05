/*This is the file to get the different individuals selected*/
void nselect(population *old_pop_ptr,population *pop2_ptr);

void nselect(population *old_pop_ptr,population *pop2_ptr)
{
  int *fit_ptr1,*fit_ptr2;
  float rnd2,*f1_ptr,*f2_ptr;
  int *s1_ptr,*s2_ptr,*select_ptr,*ptry;
  void *j,*j1;
  int i,rnd,k,n;
  
  for(n = 0,k = 0; n < popsize; n++,k++)
    {
      pop2_ptr->ind_ptr = &(pop2_ptr->ind[k]);
      select_ptr = &(pop2_ptr->ind_ptr->genes[0]);
      
      /*Select first parent randomly*/	
      rnd2 = randomperc(); 
      rnd2 = popsize* rnd2; 
      rnd = floor(rnd2);
      
      if(rnd == 0)
	rnd = popsize - k;
      if(rnd == popsize)
	rnd = (popsize-2)/2;
      
      j = &(old_pop_ptr->ind[rnd-1]);
      
      /*Select second parent randomly*/ 
      rnd2 = randomperc(); 
      rnd2 = popsize * rnd2; 
      rnd = floor(rnd2);
      
      if (rnd == 0)
	rnd = popsize - n;
      if(rnd == popsize)
	rnd = (popsize - 4)/2;
           
      j1 = &(old_pop_ptr->ind[rnd-1]);
      
      old_pop_ptr->ind_ptr = j;
      
      s1_ptr = &(old_pop_ptr->ind_ptr->genes[0]);
      fit_ptr1 = &(old_pop_ptr->ind_ptr->rank);
      f1_ptr = &(old_pop_ptr->ind_ptr->cub_len);
      
      old_pop_ptr->ind_ptr = j1;
      s2_ptr = &(old_pop_ptr->ind_ptr->genes[0]);
      fit_ptr2 = &(old_pop_ptr->ind_ptr->rank);
      f2_ptr = &(old_pop_ptr->ind_ptr->cub_len);
     //SELECTION PROCEDURE
      /*Comparing the fitnesses*/
      
      if(*fit_ptr1 > *fit_ptr2)
	{
	  for(i = 0;i < chrom;i++)
	    *select_ptr++=*s2_ptr++;
	}
      else if(*fit_ptr1 < *fit_ptr2)
	{
	  for(i = 0;i < chrom;i++)
	    *select_ptr++=*s1_ptr++;
	}
      else
	{
	  if(*f1_ptr < *f2_ptr)
	    {
	      for(i = 0;i < chrom;i++)
		*select_ptr++=*s2_ptr++;
	    }
	  else
	    {
	      for(i = 0;i < chrom;i++)
		*select_ptr++=*s1_ptr++;
	    }
	}
    }
  
return;
}








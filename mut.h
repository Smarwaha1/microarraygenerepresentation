/* This is the module used to formulate the mutation routine*/
void mutate(population *new_pop_ptr);             

void mutate(population *new_pop_ptr)
{
  int i,j,r,r1,flag,tmp[optima];
  float rand1;
  
  for(j = 0; j < popsize; j++)
    {
      new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);     
      flag = 0;
      
      for(i = 0; i < chrom; i++)
	{
	  if(new_pop_ptr->ind_ptr->genes[i] == 1)
	    {
	      tmp[flag] = i;	  
	      flag++;
	    }
	}

      for(i = 0; i < flag; i++)
	{
	  rand1 = randomperc();
	  
	  if(rand1 <= pmut_b*flag)
	    {
	      r1 = rnd(1, flag);
	      new_pop_ptr->ind_ptr->genes[tmp[r1-1]] = 0;
	      
	      do{ 
		r = rnd(1, chrom);
	      }while(new_pop_ptr->ind_ptr->genes[r-1] == 1);
	      
	      new_pop_ptr->ind_ptr->genes[r-1] = 1;
	      tmp[r1-1] = r-1;
	    }
	}
    }
  return;
}

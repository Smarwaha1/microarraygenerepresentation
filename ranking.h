/*This also demarkates the different Pareto Fronts*/
void ranking(population *pop_ptr);
int indcmp(float *ptr1,float *ptr2);

void ranking(population *pop_ptr)
{
  int i,j,k;                                 //counters
  int rnk,val,nondom,rankarr[maxpop],q;
  float *ptr1,*ptr2;
// RANKING 
  rnk = 0;  nondom = 0 ;
  
  /*Initializing all the flags to 2*/
  for(i = 0; i < popsize; i++)
    pop_ptr->ind[i].flag = 2;
  
  for(k = 0,q = 0; k < popsize; k++)
    {
      for(j = 0;j < popsize;j++)
	{
	  if(pop_ptr->ind[j].flag != 1)
	    break;
	}  /*Break if all the individuals are assigned a rank*/
      
      if(j == popsize)
	break;
      
      rnk = rnk + 1;
      
      for( j = 0 ;j < popsize; j++)
	{
	  if(pop_ptr->ind[j].flag == 0)
	    pop_ptr->ind[j].flag = 2;
	}  /*Set the flag of dominated individuals to 2*/
      
      for(i = 0;i < popsize ; i++)
	{
	  /*Select an individual which rank to be assigned*/
	  pop_ptr->ind_ptr = &(pop_ptr->ind[i]);
	  
	  if(pop_ptr->ind_ptr->flag != 1 && pop_ptr->ind_ptr->flag != 0) 
	    {
	      ptr1 = &(pop_ptr->ind_ptr->fitness[0]);
	      
	      for(j = 0;j < popsize ; j++)
		{
		  /*Select the other individual which has not got a rank*/
		  if( i != j)
		    {
		      if(pop_ptr->ind[j].flag != 1)
			{
			  pop_ptr->ind_ptr = &(pop_ptr->ind[j]);
			  ptr2 = &(pop_ptr->ind_ptr->fitness[0]);
			  
			  /*Compare the two individuals for fitness*/
			  val = indcmp(ptr1,ptr2);
			  
			  if( val == 2)
			    { 
			      pop_ptr->ind[i].flag = 0; /*individual 1 is dominated */
			      break;
			    }
			  
			  else if(val == 1)
			    pop_ptr->ind[j].flag = 0;   /*individual 2 is dominated */
			  
			  else if(val == 3)
			    {
			      nondom++;         /*individual 1 & 2 are non dominated */
			      if(pop_ptr->ind[j].flag != 0)
				pop_ptr->ind[j].flag = 3;
			    }
			  
			}  
		    }   
		}          
	      
	      if(j == popsize)
		{
		  /*Assign the rank and set the flag*/
		  pop_ptr->ind[i].rank = rnk;
		  pop_ptr->ind[i].flag = 1;
		  rankarr[q] = i;
		  q++;
		}
	    }     
	}          
      pop_ptr->rankno[rnk-1] = q;
    }
  
  pop_ptr->maxrank = rnk;
  
  return;
}


  //Routine Comparing the two individuals

int indcmp(float *ptr1,float *ptr2)
{
  float fit1[maxfun],fit2[maxfun];
  float sum23_1,sum23_2;
  int i,value,m,n;

  for(i = 0;i < nfunc ;i++)
    {
      fit1[i] = *ptr1++;
      fit2[i] = *ptr2++;
    }
  
  sum23_1 = fit1[0] + fit1[1];
  sum23_2 = fit2[0] + fit2[1];

  if(fabs(sum23_1 - sum23_2) < 1e-5) {
    if ((fit1[2] < fit2[2]) && (fit2[2]-fit1[2] <= 10))  
      { value = 3; return(value); }
    else if ((fit2[2] < fit1[2]) && (fit1[2]-fit2[2] <= 10))
      { value = 3; return(value); }
  }

 
  m = 0;
  n = 0;
  
  while(m < nfunc && fit1[m] <= fit2[m]) 
    {
      if((fit2[m] -  fit1[m]) < 1e-7) n++;
      m++;
    }
  
  if(m == nfunc) 
    {
      if(n == nfunc) value = 3;
      else value = 1;                /*value = 1 for dominationg*/
    }
  else 
    {
      m = 0;
      n = 0;
      while(m < nfunc && fit1[m] >= fit2[m])
	{
	  if((fit1[m] - fit2[m]) < 1e-7) n++;
	  m++;
	}
      
      if(m == nfunc)
	{
	  if(n != nfunc)
	    value = 2;                /*value =  2 for dominated */
	  else value =3;
	}
      else value = 3;                 /*value = 3 for incomparable*/
    }
 
  return value;
}





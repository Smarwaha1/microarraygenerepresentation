void func(population *pop_ptr);


void func(population *pop_ptr)
{
  int i, k, f[maxfun]; 
  
  pop_ptr->maxrank = 0; /*Initializing the max rank to zero*/
  
  for(i = 0; i < popsize; i++)
    {         
      pop_ptr->ind_ptr = &(pop_ptr->ind[i]);/*Population initialization*/
      
      /*calculate f1, f2, and f3 of an individual*/
      objective(pop_ptr->ind_ptr,f);  
      
      /*First objective function: mismatches in training samples*/ 
      pop_ptr->ind_ptr->fitness[0] = (float)(f[0]); // /(float)(no_train_sample);
      /*second objective function: mismatches in test samples */
      pop_ptr->ind_ptr->fitness[1] = (float)(f[1]); // /(float)(no_test_sample1);
      /*Third objective function: gene subset size*/
      pop_ptr->ind_ptr->fitness[2] = (float)(f[2]); // /(float)(chrom);
      /*End of defining objetives*/     
    }
 // RANKING
  ranking(pop_ptr);
  
  return;
}

 // Calculating the individual fitness values

objective(individual *pop_ptr_ind, int *f)
{
  int i,j,k,l,ii,pred_class,flag;
  float v1[MAXSMP],v2[MAXSMP],p,ps;
  
  f[0]=0; f[1]=0;
  flag = no_train_sample + no_test_sample1;
  
  for(ii = 0; ii < MAXSL; ii++)
    {  
      set_ptr = &(setno[ii]);    f[2]=0;
      
      for(j = 0; j < no_sample; v1[j]=0.0,v2[j]=0.0,j++); 
      
      for(i = 0; i < chrom; i++)
	{
	  if(pop_ptr_ind->genes[i] == 1)
	    {
	      set_ptr->g_ptr = &(set_ptr->gene[i]);
	      f[2] = f[2] + 1;
        		      
	      for(j = 0; j < flag; j++)
		{
		  if(j < no_train_sample)
		    { l=set_ptr->trainset[j]; k=j; }
		  else 
		    { l=set_ptr->testset1[j-no_train_sample]; k=no_train_sample; }
		  
		  p = set_ptr->g_ptr->weighting_factor[k]*
		    (set_ptr->g_ptr->exp_level[l] - set_ptr->g_ptr->weighted_vote[k]);
		  
		  if(p >= 0.0)
		    v1[l] = v1[l] + p; 
		  else
		    v2[l] = v2[l] + fabs(p);
		}
	    }
	}
          
      for(i = 0; i < flag; i++)
	{
	  if(i < no_train_sample)
	    l = set_ptr->trainset[i];	   
	  else
	    l = set_ptr->testset1[i-no_train_sample]; 
	  
	  if(v1[l] > v2[l])  //prediction strength calculation
	    {
	      ps = (v1[l]-v2[l])/(v1[l]+v2[l]);		  
	      pred_class = 1;
	    }
	  else if(v2[l] > v1[l])
	    {
	      ps = (v2[l]-v1[l])/(v1[l]+v2[l]);		  
	      pred_class = 2;
	    }
	  else
	    ps = 0.0;
	  
	  if(ps <= epsilon) //prediction strength threshold checking
	    pred_class = 0;
	  	
	  if(pred_class != actual_sample_class[l]) //Compare predicted with actual class
	    {
	      if(i < no_train_sample)
		f[0] = f[0] + 1;
	      else
		f[1] = f[1] + 1;
	    }
	}
    }
}


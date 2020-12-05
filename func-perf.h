
int fun_perf_esti(individual *ind_ptr);

int fun_perf_esti(individual *ind_ptr)
{
  int i,j,k,l,ii,pred_class,flag,count;
  float v1[MAXSMP],v2[MAXSMP],p,ps;
  
  count = 0;
  k = no_train_sample;
  flag = no_sample - no_test_sample1 - no_train_sample; 

  for(ii = 0; ii < MAXSL; ii++)
    {  
      set_ptr = &(setno[ii]);   

      for(j = 0; j < flag; v1[j]=0.0,v2[j]=0.0,j++); 
      
      for(i = 0; i < chrom; i++)
	{
	  if(ind_ptr->genes[i] == 1)
	    {
	      set_ptr->g_ptr = &(set_ptr->gene[i]);
	              		      
	      for(j = 0; j < flag; j++)
		{
		  l = set_ptr->testset2[j];  
		  
		  p = set_ptr->g_ptr->weighting_factor[k]*
		    (set_ptr->g_ptr->exp_level[l] - set_ptr->g_ptr->weighted_vote[k]);
		  
		  if(p >= 0.0)
		    v1[j] = v1[j] + p; 
		  else
		    v2[j] = v2[j] + fabs(p);
		}
	    }
	}
          
      for(i = 0; i < flag; i++)
	{
	  l = set_ptr->testset2[i];

	  if(v1[i] > v2[i])  //prediction strength calculation
	    {
	      ps = (v1[i]-v2[i])/(v1[i]+v2[i]);		  
	      pred_class = 1;
	    }
	  else if(v2[i] > v1[i])
	    {
	      ps = (v2[i]-v1[i])/(v1[i]+v2[i]);		  
	      pred_class = 2;
	    }
	  else
	    ps = 0.0;
	  
	  if(ps <= epsilon) 
	    pred_class = 0;
	  	
	  printf("l=%d,pred=%d,actual=%d\n",l,pred_class,actual_sample_class[l]); 	  
	 
	  if(pred_class != actual_sample_class[l]) 
	    count++;	      
	}
    }
  return count;
}
/*===========================================================================*/

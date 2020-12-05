
void keepalive(population *pop1_ptr,population *pop2_ptr,population *pop3_ptr,int gen);

typedef struct
{
  int maxrank,                    /*Max rank of the global population*/
    rankar[2*maxpop][2*maxpop],   /*array of individual numbers at a particular rank */
    rankno[2*maxpop];             /*record of no. of individuals at a particular rank*/
  int genes[2*maxpop][maxchrom],      
    rank[2*maxpop],               /*rank of different individuals*/
    flag[2*maxpop];               /*Setting the flag */
  float fitness[2*maxpop][maxfun],/*Fitness values for the different individuals*/
    cub_len[2*maxpop];            /*Dummyfitness*/
}globpop;
/*Population structure for the pool having both the old as well as new population*/
globpop globalpop,*global_pop_ptr;

/*Ranking the global pool*/
void grank(int gen);
/*Comparison of the variables*/
int indcmp1(float *ptr1,float *ptr2);
/*Sorting for the function values in ascending order*/
void gsort(int rnk,int sel);
/* Fitness crowding distance calculation */
void fitness_crowding(int rnk);
/*Multi-Modal Fitness crowding*/
void crowding_fitness(int rnk, int sel);
/*Duplicates elimination for Multi-Modal*/
void duplicate_elimination(int rnk);
/*Sorting the function values in ascending order*/
void sort(int rnk);
/*Initialization of the global variables*/
int left,Lastrank;
float fpara1[2*maxpop][2];

/*Function starts here*/
void keepalive(population *pop1_ptr,population *pop2_ptr,population *pop3_ptr,int gen)
{
  int i,j,jj,k,m,a1,l,front_pop[maxpop];
  int sum,st,str,pool,poolf,sel,r1;
  int *gene1_ptr, *gene2_ptr,leftsum;
  float rnd,a,*gene3_ptr,x3,*gene4_ptr;
  /*----------------------------* RANKING *---------------------------------*/
  
  for(i = 0; i < popsize; i++)
    {
      /*Binary Coded GA genes are copied*/
      for(k = 0; k < chrom; k++)
	{
	  globalpop.genes[i][k]=pop1_ptr->ind[i].genes[k];
	  globalpop.genes[i+popsize][k] = pop2_ptr->ind[i].genes[k];
	}
      /*Fitness is copied to the global pool*/
      for(l = 0; l < nfunc; l++)
	{
	  globalpop.fitness[i][l] = pop1_ptr->ind[i].fitness[l];
	  globalpop.fitness[i+popsize][l] = pop2_ptr->ind[i].fitness[l];
	}
      /*Initialising the dummy fitness to zero*/
      globalpop.cub_len[i] = 0;
      globalpop.cub_len[i+popsize] = 0;
    }
    
  /*Finding the global ranks */
  global_pop_ptr = &(globalpop);  
  grank(gen);
  
  m = globalpop.maxrank;
  pool = 0;
  
  /*Initializing the flags of population to zero */
  for(i = 0; i < 2*popsize; i++)
    globalpop.flag[i] = 0;
  
  /*Fitness Crowding*/
  for(i = 0; i < m; i++)
    fitness_crowding(i+1);
  
  /*Finding number of solutions are allowed to copy for pop3
   from each front by eliminating the Duplicates for multimodal */
  for(i = 0; i < m; i++)
    {
      duplicate_elimination(i+1); /*Duplicate Elimination*/
      st = pool;
      pool+= globalpop.rankno[i];
      
      if(pool <= popsize)
	{
	  for(k = 0; k < 2*popsize; k++)
	    if(globalpop.rank[k] == i+1)
	      globalpop.flag[k] = 1;
	  
	  pop3_ptr->rankno[i] = globalpop.rankno[i];
	  
	  if(pool == popsize) break;
	}
      else
  	{	
	  sel = popsize - st;
	  Lastrank = i+1;
	  pop3_ptr->rankno[i] = sel;
	  crowding_fitness(i+1, sel); //Multi-Modal NSGA-II
	  break;
	}
    }
  
 
  
  for(i = 0,k = 0;i < 2*popsize && k < popsize; i++)
    {
      if(globalpop.flag[i] == 1)
	{
	  gene1_ptr = &(globalpop.genes[i][0]);
	  pop3_ptr->ind_ptr = &(pop3_ptr->ind[k]);
	  gene2_ptr = &(pop3_ptr->ind_ptr->genes[0]);
	  
	  for(j = 0 ; j < chrom; j++)
	    *gene2_ptr++ = *gene1_ptr++;
	  for(j = 0;j < nfunc;j++)
	    pop3_ptr->ind[k].fitness[j] = globalpop.fitness[i][j];
	  pop3_ptr->ind[k].cub_len = globalpop.cub_len[i];
	  pop3_ptr->ind[k].rank = globalpop.rank[i];
	  k++;  // increment the pop3 counter
	}
    }
  pop3_ptr->maxrank = Lastrank;
  return;
}

//Ranking the Fitness values in ascending order

void grank(int gen)
{
  int i,j,k,rnk,val,nondom,popsize1,gflg[2*maxpop],q;
  float *ptr1,*ptr2;
  
  rnk = 0;
  nondom = 0;
  popsize1 = 2*popsize;
  
  for(i = 0;i < popsize1;i++)
    gflg[i] = 2;
  
  for(k = 0;k < popsize1;k++)
    {
      q =  0;
      for(j = 0;j < popsize1;j++)
	if (gflg[j] != 1) break;
      
      if(j == popsize1) break;
      rnk = rnk +1;
      for( j = 0 ;j < popsize1; j++)
	if(gflg[j] == 0) gflg[j] = 2;
      
      for(i = 0;i < popsize1 ; i++)
	{
	  if(gflg[i] != 1 && gflg[i] != 0) 
	    {
	      ptr1 = &(global_pop_ptr->fitness[i][0]);
	      for(j = 0;j < popsize1 ; j++)
		{
		  if( i!= j)
		    {
		      if(gflg[j] != 1)
			{
			  ptr2 = &(global_pop_ptr->fitness[j][0]);
			  val = indcmp1(ptr1,ptr2);
			  if( val == 2)
			    { 
			      gflg[i] = 0;/* individual 1 is dominated */
			      break;
			    }
			  if(val == 1)
			    {
			      gflg[j] = 0;/* individual 2 is dominated */
			    }
			  if(val == 3)
			    {
			      nondom++;/* individual 1 & 2 are non dominated */
			      if(gflg[j] != 0)gflg[j] = 3;
			    }
			}
		    }
		}
	      if( j == popsize1)
		{
		  global_pop_ptr->rank[i] = rnk;
		  gflg[i] = 1;
		  global_pop_ptr->rankar[rnk-1][q] =  i;
		  q++;
		}
	    }
	}
      global_pop_ptr->rankno[rnk-1] = q;
    } 
  global_pop_ptr->maxrank = rnk;
  return;
}

//Two Individuals are compared by non dominated principles

int indcmp1(float *ptr1,float *ptr2)
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
  
  m = 0; n = 0;
  while(m < nfunc && fit1[m] <= fit2[m])
    {
      if((fit2[m] - fit1[m]) < 1e-7) n++;
      m++;
    }
 
  if(m == nfunc) 
    {
      if(n == nfunc) value = 3; /*value = 3 for incomparable*/
      else value = 1; /*value = 3 for incomparable*/
    }
  else 
    {
      m = 0;n = 0;
      while(m < nfunc && fit1[m] >= fit2[m]) 
	{
	  if((fit1[m] - fit2[m]) < 1e-7) n++;
	  m++;
	}
      
      if(m == nfunc)
	{
	  if(n != nfunc) value = 2;   /*value =  3 for incomparable */
	  else value = 3;  /*value = 3 for incomparable*/
	}
      else value = 3;                   /*value = 3 for incomparable*/
    }
  return value;
}

 //This is the file used to sort the dummyfitness arrays

void gsort(int rnk,int sel)
{
  int i,j,a,q;
  float array[2*maxpop][2],temp,temp1;
  
  q = globalpop.rankno[rnk-1];
  
  for(i = 0 ;i < q ;i++)
    {
      array[i][0] = globalpop.rankar[rnk-1][i];
      a = globalpop.rankar[rnk-1][i];
      array[i][1] = globalpop.cub_len[a];
    }
  for(i = 0;i < q ;i++)
    {
      for(j = i+1;j < q;j++)
	{
	  if(array[i][1] < array[j][1])
	    {
	      temp = array[i][1];
	      temp1 = array[i][0];
	      array[i][1] = array[j][1];
	      array[i][0] = array[j][0];
	      array[j][1] = temp;
	      array[j][0] = temp1;
	    }
	}
    }
  
  for(i = 0;i < sel;i++)
    {
      a = array[i][0];
      globalpop.flag[a] = 1;
    }
  return;
}

//Assigning dummy fitness values to the Global population members

void fitness_crowding(int rnk)
{
  float length[2*maxpop][2],max;
  int i,j,m1,a ;

  m1 = globalpop.rankno[rnk-1];
  
  for(j = 0;j < nfunc;j++)
    {
      for(i = 0;i < m1;i++)
	{
	  fpara1[i][0] = 0;
	  fpara1[i][1] = 0;
	}
      
      for(i = 0;i < m1;i++)
	{
	  a = globalpop.rankar[rnk-1][i];
	  fpara1[i][0] = (float)a ;
	  fpara1[i][1] = globalpop.fitness[a][j];
	}
      
      sort(m1); /*Sort the arrays in ascending order of the fitness*/
      
      max = fpara1[m1-1][1];
      for(i = 0;i < m1;i++)
	{
	  if(i == 0 ||i == (m1-1))
	    { 
	      length[i][0] = fpara1[i][0];
	      length[i][1] = 100*max;
	    }
	  else
	    {
	      length[i][0] = fpara1[i][0];
	      length[i][1] = fabs(fpara1[i+1][1]- fpara1[i-1][1]);
	    }
	}
      for(i = 0;i < m1;i++)
	{
	  a = length[i][0];
	  globalpop.cub_len[a] += length[i][1];
	}
    }
  return;
}

  //Sorting the fitness values in ascending order

void sort(int m1)
{
  float temp,temp1; 
  int i1,j1,k1;
  for(k1 = 0;k1 < m1-1;k1++)
    {
      for(i1 = k1+1;i1 < m1;i1++)
	{
	  if(fpara1[k1][1] > fpara1[i1][1])
	    {
	      temp = fpara1[k1][1];
	      temp1 = fpara1[k1][0];
	      fpara1[k1][1] = fpara1[i1][1];
	      fpara1[k1][0] = fpara1[i1][0];
	      fpara1[i1][1] = temp;
	      fpara1[i1][0] = temp1;
	    }
	}
    }
  return;
}

 //Subroutine to eliminate the duplicates 

void duplicate_elimination(int rnk)
{
  int i,j,k,m1,a1,a2,q,m;
  int *gene1_ptr,*gene2_ptr,*flag;
  
  m1 = globalpop.rankno[rnk-1];
  q=0;

  flag = (int *)calloc(m1, sizeof(int));
 
  for(i = 0; i < (m1-1); i++)
    {
      if(flag[i] == 0)
	{
	  a1 = globalpop.rankar[rnk-1][i];
	  for(j = i+1; j < m1; j++)
	    {
	      if(flag[j]==0)
	      	{
		  gene1_ptr = &(globalpop.genes[a1][0]);
		  a2 = globalpop.rankar[rnk-1][j];  
		  gene2_ptr = &(globalpop.genes[a2][0]);
		  m = 0;
		  while((m < nfunc) && (globalpop.fitness[a1][m] == globalpop.fitness[a2][m])) 
		    {
		      m++;
		    }
		  
		  if(m == nfunc)
		    {
		      k=0;
		      while((k < chrom) && (*gene1_ptr++ == *gene2_ptr++))
			{k++;}
		      if(k == chrom)
			{
			  globalpop.rank[a2] = globalpop.maxrank + 1;
			  flag[j] = 1;
			  q++;
			}
		    }
		}
	    }
	}
    }
  
  for(i = 0,k = 0; i < m1; i++)
    {
      if(flag[i] == 0)
	{
	  globalpop.rankar[rnk-1][k] = globalpop.rankar[rnk-1][i];
	  k++;
	}
    }
  
  globalpop.rankno[rnk-1] = globalpop.rankno[rnk-1] - q;
  free(flag);
  return;
}

//This is the file used to sort the dummyfitness arrays

void crowding_fitness(int rnk, int sel)
{
  int i,j,k,a1,a2,q,m,*flag;
  int dist_arr[2*maxpop][maxpop],dist_no[2*maxpop],dist_no1[2*maxpop],count;
  float a;

  q = globalpop.rankno[rnk-1];
  flag = (int *)calloc(q, sizeof(int));
  
  for(i = 0, count = 0; i < q; i++)
    {
      a1 = globalpop.rankar[rnk-1][i];
      globalpop.flag[a1] = 0;
      
      if(flag[i] == 0)
	{
	  k = 0;
	  dist_no[count] = 1;
	  dist_arr[count][k] = a1;
	  
	  for(j = i+1; j < q; j++)
	    {
	      if(flag[j] == 0)
		{
		  a2 = globalpop.rankar[rnk-1][j];
		  m = 0;
		  while((m < nfunc) && (globalpop.fitness[a1][m] == globalpop.fitness[a2][m])) 
		    {
		      m++;
		    }
		  
		  if(m == nfunc)
		    {
		      k = k+1;
		      flag[j] = 1;
		      dist_no[count] = dist_no[count]+1;
		      dist_arr[count][k] = a2;
		    }
		}
	    }
	  count = count + 1;
	}
    }

 free(flag);

 if(count >= sel)
    gsort(rnk,sel);      
 else
   {
     m = 0;
     for(i = 0; i < count; i++)
       {
	 a = (float)((dist_no[i]-1)*(sel - count))/(float)(q - count);
	 dist_no1[i] = a + 1.4999;
	 m = m + dist_no1[i];
       }
     do
       {
	 if(m < sel)
	   {
	     a1 = rnd(1, count);
	     if( dist_no[a1-1] > dist_no1[a1-1])
	       {
		 dist_no1[a1-1] = dist_no1[a1-1]+1;
		 m = m + 1;
	       }
	   }
	 else if(m > sel)
	   {
	     a1 = rnd(1, count);
	     if(dist_no1[a1-1] > 1)
	       {
		 dist_no1[a1-1] = dist_no1[a1-1]-1;
		 m = m - 1;
	       }
	   }
	 
	 if(m == sel)
	   {
	     for(i = 0; i < count; i++)
	       {
		 for(j = 0; j < dist_no1[i]; j++)
		   {
		     a1 = dist_arr[i][j];
		     globalpop.flag[a1] = 1;
		   } 
	       }
	   }
       } while(m != sel);  
   }
 return;
}


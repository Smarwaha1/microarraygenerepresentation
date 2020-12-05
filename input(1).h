/*This is a file to get the input for the GA program*/
void input(FILE *rep_ptr);

void input(FILE *rep_ptr)
{
  /*Gene Expression Dataset Specifications*/ 
  printf("\nSpecify Gene expression dataset Parameters (Please be carefull)->\n");

  /*Asks for number of genes present in the dataset*/
  printf("\nGive the number of genes (see the dataset) :");
  scanf("%d", &chrom);
  
 
  if(chrom > maxchrom)
    {
      printf("Increase maxchrom in nsga2.c, Currently set :%d\n",maxchrom);
      exit(1);
    }
  
  /*Asks for number of samples present in the dataset*/
  printf("\nGive the number of samples (see the dataset) :");
  scanf("%d", &no_sample);
  
  
  if(no_sample > MAXSMP)
    {
      printf("Increase MAXSMP in nsga2.c, Currently set :%d\n",MAXSMP);
      exit(1);
    }
  
  /*Asks for number of training samples used to construct the classifier*/
  printf("\nGive the number of training samples (>50 perc of total samples) :");
  scanf("%d", &no_train_sample);
  
   
  if(no_train_sample > MAXTR)
    {
      printf("Increase MAXTR in nsga2.c, Currently set %d.\n",MAXTR);
      exit(1);
    }

  /*NSGA-II Parameters*/
  printf("\nSpecify NSGA-II Parameters (Please be carefull)->\n");
  
  
  printf("\nGive the no.of generations (100 to 1000) :");
  scanf("%d",&gener);
  
  /*Asks for number of the individuals in the population*/
  printf("\nGive the Population size (an even no. 100 to 1000) :");
  scanf("%d",&popsize);


  if(popsize > maxpop)
    {
      printf("Increase maxpop in nsga2.c, Currently set %d\n",maxpop);
      exit(1);
    }
    
  /*Asks for number of the objective functions*/
  printf("\nGive no. of objective functions (fixed 3) :");
  scanf("%d",&nfunc);


  if(nfunc > maxfun)
    {
      printf("Sorry! Here you are recommanded to use 3 for this classification problem\n");
      exit(1);
    }

  /*Asks for the Crossover probability*/
  printf("\nGive the cross-over probability (between 0.5 and 1) :");
  scanf("%f",&pcross);

  /*Asks for the Mutation probability*/
  printf("\nGive the mutation probability for binary strings (between 0 and 1/l) :");
  scanf("%f",&pmut_b);

  /*Asks for the random seed for generating random numbers*/   
  printf("\nGive random seed (between 0 and 1) :");
  scanf("%f",&seed);
 
  /*Print the GA parameters in the file output.dat */
  fprintf(rep_ptr,"==================================================\n");
  fprintf(rep_ptr,"\tGene Expression Dataset Specifications\n");
  fprintf(rep_ptr,"==================================================\n");
  fprintf(rep_ptr,"No. of genes in dataset :%d\n",chrom);
  fprintf(rep_ptr,"No. of samples in dataset :%d\n",no_sample);
  fprintf(rep_ptr,"No. of training samples :%d\n",no_train_sample);

  fprintf(rep_ptr,"\n==================================================\n");
  fprintf(rep_ptr,"\tNSGA-II Parameters\n");
  fprintf(rep_ptr,"==================================================\n");
  fprintf(rep_ptr,"No. of generations ->%d\n",gener); 
  fprintf(rep_ptr,"Population Size ->%d\n",popsize);
  fprintf(rep_ptr,"Chromosome Length ->%d\n",chrom);
  fprintf(rep_ptr,"No. of Functions ->%d\n",nfunc);
  fprintf(rep_ptr,"X-over on binary string is SINGLE POINT X-OVER\n");
  fprintf(rep_ptr,"Cross-over Probability ->%f\n",pcross);
  fprintf(rep_ptr,"Mutation Probability for binary strings -> %f\n",pmut_b);
  fprintf(rep_ptr,"Random Seed ->%f\n",seed);
  
  return;
}





















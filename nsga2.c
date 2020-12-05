
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define square(x) ((x)*(x))
#define maxpop   1000       /*Max population*/
#define maxchrom 5000       /*Max chromosome length*/
#define maxfun   3          /*Max no. of functions*/
#define epsilon  0.00       /*Thresold on prediction strength*/
#define MAXSMP   100        /*Max number of available samples*/
#define MAXTR    80         /*Max number of training samples*/
#define MAXSL    1          /*Max Multiple training sets 'H'*/
#define optima   5000	    /*Optimal Gene subset size */
#define no_test_sample1  14

int gener,                /*No of generations*/
  nfunc,                  /*No of functions*/
  popsize,                /*Population Size*/
  chrom,                  /*Chromosome size*/
  sharespace;             /*Sharing space (either parameter or fitness)*/

float seed,               /*Random Seed*/
  pcross,                 /*Cross-over Probability*/
  pmut_b,                 /*Mutation Probability*/
  delta_fit,              /*Variables required for fitness for fitness sharing */
  min_fit,
  front_ratio;

time_t starttime,endtime; /*Time count variables*/ 

typedef struct            /*individual properties*/
{
  int genes[maxchrom],    /*bianry chromosome*/
    rank,                 /*Rank of the individual*/
    flag;                 /*Flag for ranking*/
  float fitness[maxfun],  /*Fitness values */
    cub_len;              /*crowding distance of the individual*/
}individual;              /*Structure defining individual*/

typedef struct            /*Population properties in a generation*/
{
  int maxrank;            /*Maximum rank present in the population*/
  float rankrat[maxpop];  /*Rank Ratio*/
  int rankno[maxpop];     /*Individual at different ranks*/
  individual ind[maxpop], /*Different Individuals*/
    *ind_ptr; 
}population;              /*Popuation Structure*/

/*Defining the population Structures*/
population oldpop,
  newpop,
  matepop,
  *old_pop_ptr,
  *new_pop_ptr,
  *mate_pop_ptr;

/*Including the necessary header files*/ 
#include "read1.h"        /*Perform preprocessing steps on expression values*/
#include "input.h"        /*File Takes Input to NSGA-II from user*/
#include "random.h"       /*Random Number Generator*/
#include "init.h"         /*Random Initialization of the population*/
#include "ranking.h"      /*File creating the pareto fronts */ 
#include "func-con.h"     /*File Having the Objective Functions*/
#include "func-perf.h" 
#include "select.h"       /*File for Crowded Tournament Selection*/
#include "crossover.h"    /*Binary Uniform Cross-over*/
#include "mut.h"          /*Binary Bit-wise Mutation*/
#include "keepaliven.h"   /*File For Elitism and Sharing Scheme*/

//Main program

main()
{
  /*Some Local variables to this Problem (Counters And some other pointers*/
  int i,j,k,l,f,maxrank1,test_error;
  float *ptr,tot;
  /*File Pointers*/
  FILE *rep_ptr, *end_ptr, *g_var;
  
  /*Initialization of file pointers*/  
  rep_ptr = fopen("output.inp","w");     
  end_ptr = fopen("final_fitness.inp","w");
  g_var = fopen("final_var.inp","w");      
  
  /*Get the input from the file input.h*/
  input(rep_ptr);

  /*Initialize the random no generator*/
  warmup_random(seed);

  /*Calling functions to read expression levels*/
  read_dataset();
  calc_mean_sd();   
  time(&starttime);

  /*Population initialization*/
  old_pop_ptr = &(oldpop); 
  init(old_pop_ptr);  
  
  /*Population Assignment*/
  old_pop_ptr = &(oldpop);
  new_pop_ptr = &(newpop);
  
  for(j = 0;j < popsize;j++)
    {
      old_pop_ptr->rankno[j] = 0;
      old_pop_ptr->ind[j].cub_len = 0;
      new_pop_ptr->rankno[j] = 0;
      new_pop_ptr->ind[j].cub_len = 0;
    }

  /*Function called for fitness calculation*/
  old_pop_ptr = &(oldpop);     
  func(old_pop_ptr);

  /*Function called for fitness calculation*/
  old_pop_ptr = &(oldpop); 
  ranking(old_pop_ptr);
  
  /*GENERATION STARTS HERE*/
  for(i = 0; i < gener; i++)
    {         
      printf("Generation No: %d\n",i+1);

      /*SELECTION*/
      old_pop_ptr = &(oldpop);
      mate_pop_ptr = &(matepop);
      nselect(old_pop_ptr, mate_pop_ptr);
	  
      /*BINARY CROSSOVER*/      
      new_pop_ptr = &(newpop);
      mate_pop_ptr = &(matepop);
      crossover(new_pop_ptr, mate_pop_ptr );
      
      /*BINARY MUTATION*/
      new_pop_ptr = &(newpop);
      mutate(new_pop_ptr);

      /*FUNCTION EVALUATION*/
      new_pop_ptr = &(newpop);
      func(new_pop_ptr);
      
      /*SELECTION KEEPING FRONTS ALIVE*/
      old_pop_ptr = &(oldpop);
      new_pop_ptr = &(newpop);
      mate_pop_ptr = &(matepop);
	  
      /*Elitism And Sharing Implemented*/
      keepalive(old_pop_ptr ,new_pop_ptr ,mate_pop_ptr,i+1);      

      /*Rank Ratio Calculation*/
      mate_pop_ptr = &(matepop);
      new_pop_ptr = &(matepop);
      old_pop_ptr = &(oldpop);
      
      /*Finding the greater maxrank among the two populations*/
      if(old_pop_ptr->maxrank > new_pop_ptr->maxrank)
	maxrank1 = old_pop_ptr->maxrank;
      else 
	maxrank1 = new_pop_ptr->maxrank;

      for(j = 0;j < maxrank1 ; j++)
	{ 
	  /*Sum of the no of individuals at any rank in old and new populaion*/
	  tot = (old_pop_ptr->rankno[j])+ (new_pop_ptr->rankno[j]);
	  /*Finding the rank ratio for new population at this rank*/
	  new_pop_ptr->rankrat[j] = (new_pop_ptr->rankno[j])/tot;
	}
      
      //Copying the new population to old population
      old_pop_ptr = &(oldpop);
      new_pop_ptr = &(matepop);
      
      for(j = 0; j < popsize; j++)
	{
	  old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[j]);
	  new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[j]);
	  
	  /*For Binary GA copying of the chromosome*/
	  for(l = 0; l < chrom; l++)
	    old_pop_ptr->ind_ptr->genes[l] = new_pop_ptr->ind_ptr->genes[l];
	  
	  /*Copying the fitness vector */	  
	  for(l = 0; l < nfunc; l++)
	    old_pop_ptr->ind_ptr->fitness[l] = new_pop_ptr->ind_ptr->fitness[l];
	  
	  /*Copying the dummy fitness*/
	  old_pop_ptr->ind_ptr->cub_len = new_pop_ptr->ind_ptr->cub_len;
	  
	  /*Copying the rank of the individuals*/
	  old_pop_ptr->ind_ptr->rank = new_pop_ptr->ind_ptr->rank;

	  /*Copying the flag of the individuals*/
	  old_pop_ptr->ind_ptr->flag = new_pop_ptr->ind_ptr->flag;
	}//end of j

      maxrank1 = new_pop_ptr->maxrank;

      /*Copying the array having the record of the individual at different ranks */
      for(l = 0; l < popsize; l++)
	old_pop_ptr->rankno[l] = new_pop_ptr->rankno[l];

      /*Copying the maxrank */
      old_pop_ptr->maxrank = new_pop_ptr->maxrank;

    }/*end of generation loop i*/

  /*Find the feasible solutions at the end of generation and Print the results*/
  old_pop_ptr = &(matepop);
  time(&endtime);
  /*report the CPU time*/
  fprintf(rep_ptr,"\nTime taken = %lf\n",difftime(endtime, starttime));

  for(f = 0; f < popsize ; f++)
    {
      old_pop_ptr->ind_ptr = &(old_pop_ptr->ind[f]);

      if (old_pop_ptr->ind_ptr->rank == 1)
	{
	  test_error = fun_perf_esti(old_pop_ptr->ind_ptr);

	  fprintf(end_ptr,"%2.1f\t",old_pop_ptr->ind_ptr->fitness[2]);//*chrom);
	  fprintf(end_ptr,"%2.1f\t",old_pop_ptr->ind_ptr->fitness[0]);//*no_train_sample);
	  fprintf(end_ptr,"%2.1f\t",old_pop_ptr->ind_ptr->fitness[1]);//*no_test_sample1);
	  fprintf(end_ptr,"%d\t",test_error);

	  for(l = 0; l < chrom; l++)
	    {
	      if(old_pop_ptr->ind_ptr->genes[l] == 1)
		fprintf(end_ptr,"%d\t",l+1);
	    }
	  fprintf(end_ptr,"\n");
	}//feasibility check
    }//end of f (printing)

  printf("NOW YOU CAN LOOK IN THE FILE OUTPUT2.DAT\n");
  /*Closing the files*/
  fclose(rep_ptr);
  fclose(end_ptr);
  fclose(g_var);
}





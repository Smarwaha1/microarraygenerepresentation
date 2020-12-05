/*File reads the data from the input file and construct the classifier*/
/*Structure defines the properties of a gene*/
typedef struct
{  
  float exp_level[MAXSMP];       /*Expression levels*/
  float weighting_factor[MAXTR]; /*Weighting factors by LOOCV*/
  float weighted_vote[MAXTR];    /*Weighted Vote by LOOCV*/
}DATASET;
/*Structure defines the properties of dataset*/
typedef struct
{
  int trainset[MAXTR],      /*array stores the sample numbers belonging to trainset*/
    testset1[MAXTR],
    testset2[MAXTR],
    testset[MAXTR];         /*array stores the sample numbers belonging to test set*/

  DATASET gene[maxchrom],   /*Defining the genes*/
    *g_ptr;                 /*Pointer for gene structure*/ 
}SETNUMBER;
/*Defining the structure variables*/
SETNUMBER setno[MAXSL],
  *set_ptr;
/*Initializing the Global variables*/
int no_sample,                 /*number of samples in dataset*/
  no_test_sample,              /*number of test samples in the dataset*/
  no_train_sample,             /*number of train samples*/ 
  actual_sample_class[MAXSMP]; /*array stores the class of the sample*/
int testset1[50];

void read_dataset();
void select_trainset();  
void select_testset();
void select_testset1();
void select_testset2();
void log_exp_levels(float **e_level);
void normalize_exp_levels(float **e_level); 
void calc_mean_sd();
// Function to read the data from the dataset
void read_dataset()
{
  int i, j;
  float **e_level;
  FILE *in;
  
  in = fopen("3859-gene-leukemia.txt", "rb"); // Dataset is opened
  
  if(in == NULL) 
    { printf("\nFILE IS EMPTY! cannot open the file 'input.txt'\n");
    exit(0); }
  
  for(i = 0; i < no_sample; i++)
    fscanf(in,"%d",&actual_sample_class[i]);
  
  e_level = (float **)malloc(chrom*sizeof(float*));  
  
  for(i = 0; i < chrom; i++)
    {
      *(e_level+i) = (float *)malloc(no_sample*sizeof(float));
      
      for(j = 0; j < no_sample; j++)
	fscanf(in,"%f", &e_level[i][j]);
    }
  
  
  select_trainset();  
  
select_testset();
select_testset1();
select_testset2();
log_exp_levels(e_level);   
normalize_exp_levels(e_level);  
  
  free(e_level);
  fclose(in);
  return;
}
//Function to select the train dataset 
void select_trainset()
{
  int i, ii, j;
  
  no_test_sample = no_sample - no_train_sample;
  
  for(ii = 0; ii < MAXSL; ii++)
    {      
      set_ptr = &(setno[ii]);
      
      for(i = 0; i < no_train_sample; i++)
	{
	  set_ptr->trainset[i] = i;
	  set_ptr->trainset[i] = rnd(0, (no_sample-1));
	  
	  for(j = 0; j < i; j++)
	  {
	  if (set_ptr->trainset[i] == set_ptr->trainset[j])
	  {
	  set_ptr->trainset[i] = rnd(0, (no_sample-1));
	  j = -1;
	  }
	  }
	} 
    }
  return;
}
//Function to select the test dataset
void select_testset()
{
  int i, j, ii, flag;
  
  for(ii = 0; ii < MAXSL; ii++)
    {         
      set_ptr = &(setno[ii]);  flag = 0;

      for(i = 0; i < no_sample; i++)
	{
	  for(j = 0; ((j < no_train_sample) && (i != set_ptr->trainset[j])); j++); 
	  	  
	  if(j == no_train_sample)
	    {
	      set_ptr->testset[flag] = i;
	      flag++;
	    }
	}
    }
  return;
}

//Function to select the test dataset 1
void select_testset1()
{
  int i, ii, j, debread1;
  FILE *debtest1;
  
  debtest1 = fopen("test_file1.inp","r");
  for(ii = 0; ii < MAXSL; ii++)
    {      
      set_ptr = &(setno[ii]);
      
      for(i = 0; i < no_test_sample1; i++)
	{
	   testset1[i] = rnd(0, (no_test_sample-1));
	  fscanf(debtest1,"%d",&debread1);
	   testset1[i] = debread1;
	  
	  set_ptr->testset1[i] = debread1; 
	  set_ptr->testset[testset1[i]] ;
	} 
    }
  fclose(debtest1);
return;
}

//Function to select the test dataset 2
void select_testset2()
{
  int i, j, ii, debread2;
  FILE *debtest2;
  
  debtest2 = fopen("test_file2.inp","r");
  for(ii = 0; ii < MAXSL; ii++)
    {         
      set_ptr = &(setno[ii]);
      
      for(i = 0; i < no_test_sample; i++)
	{
	  fscanf(debtest2,"%d",&debread2);
	  set_ptr->testset2[i] = debread2;
	}
    }
  fclose(debtest2);
  return;
}

//Function to claculate the logarithm of the dataset
void log_exp_levels(float **e_level)
{
  int i, j;
  
  for(i = 0; i < chrom; i++)
    {
      for(j = 0; j < no_sample; j++)
	{
	  if(e_level[i][j] < 20.0)  //flooring with 20
	    e_level[i][j] = 20.0;
	  
	  else if(e_level[i][j] > 16000.0) //ceiling with 16000
	    e_level[i][j] = 16000.0;
	  
	  e_level[i][j] = log10(e_level[i][j]); //log-transofrmation
	}
    }
  return;                                        
}
// Function to normalize the dataset
void normalize_exp_levels(float **e_level)
{
  int i, j, ii;
  float mean, sd;
  
  for(ii = 0; ii < MAXSL; ii++)
    {
      set_ptr = &(setno[ii]);

      for(i = 0; i < chrom; i++)
	{
	  mean = 0.0, sd = 0.0;      
	  set_ptr->g_ptr = &(set_ptr->gene[i]);
	  
	  for(j = 0; j < no_train_sample; j++)
	    mean = mean + e_level[i][set_ptr->trainset[j]];
	  
	  mean = mean/no_train_sample;
	  
	  for(j = 0; j < no_train_sample; j++)
	    sd = sd + pow((e_level[i][set_ptr->trainset[j]] - mean), 2);
	  
	  sd = sqrt(sd/no_train_sample);
	  
	  for(j = 0; j < no_sample; j++)
	    {
	      if(sd > 0.01)
		set_ptr->g_ptr->exp_level[j] = (e_level[i][j] - mean)/sd;
	      else
		set_ptr->g_ptr->exp_level[j] = 0.0;
	    }
	}
    }
  return;
}
//Function to calculate mean and standard deviation
void calc_mean_sd()
{
  int i, j, k, ii, flag1, flag2;
  float mean1, mean2, sd1, sd2;

  for(ii = 0; ii < MAXSL; ii++)
    {            
      set_ptr = &(setno[ii]);
      for(i = 0; i < chrom; i++)
	{
	  set_ptr->g_ptr = &(set_ptr->gene[i]);

	  for(k = 0; k <= no_train_sample; k++)
	    {
	      mean1 = 0.0;      mean2 = 0.0;
	      sd1 = 0.0;        sd2 = 0.0;
	      flag1 = 0;        flag2 = 0;
	      
	      for(j = 0; j < no_train_sample; j++)
		{
		  if(j != k)
		    {
		      if(actual_sample_class[set_ptr->trainset[j]] == 1)
			{
			  mean1 = mean1 + set_ptr->g_ptr->exp_level[set_ptr->trainset[j]];
			  flag1++;
			}
		      else
			{
			  mean2 = mean2 + set_ptr->g_ptr->exp_level[set_ptr->trainset[j]];
			  flag2++;
			}
		    }
		}
	      
	      mean1 = mean1/flag1;
	      mean2 = mean2/flag2;
	      
	      for(j = 0; j < no_train_sample; j++)
		{ 
		  if(j != k)
		    {		  
		      if(actual_sample_class[set_ptr->trainset[j]] == 1)
			sd1 = sd1 + pow((set_ptr->g_ptr->exp_level[set_ptr->trainset[j]]-mean1), 2);
		      else
			sd2 = sd2 + pow((set_ptr->g_ptr->exp_level[set_ptr->trainset[j]]-mean2), 2);
		    }
		}
	      
	      sd1 = sqrt(sd1/flag1);
	      sd2 = sqrt(sd2/flag2);
	      
	      if((sd1 == 0.0) && (sd2 == 0.0))
		set_ptr->g_ptr->weighting_factor[k] = 0.0;
	      else
		set_ptr->g_ptr->weighting_factor[k] = (mean1 - mean2)/(sd1+sd2);
	      
	      set_ptr->g_ptr->weighted_vote[k] = (mean1 + mean2)/2;
	    }
	}
    }
  return;
}


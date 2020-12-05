
static double oldrand[55];       /*Array of 55*/
static int jrand;                /*random numbers current random number*/ 
static double rndx2;             /*used with random normal deviate*/ 
static int rndcalcflag;          /*used with random normal deviate*/ 

//Create a batch 55 random numbers 
advance_random() 
{ 
  int j1; 
  double new_random; 
  
  for(j1 = 0; j1 < 24; j1++) 
    { 
      new_random = oldrand[j1] - oldrand[j1+31]; 
      
      if(new_random < 0.0) 
	new_random = new_random + 1.0; 
      
      oldrand[j1] = new_random; 
    } 
  
  for(j1 = 24; j1 < 55; j1++) 
    { 
      new_random = oldrand [j1] - oldrand [j1-24]; 
      
      if(new_random < 0.0) 
	new_random = new_random + 1.0; 
      
      oldrand[j1] = new_random; 
    } 
} 

 //Flip a biased coin - true if heads 
  
int flip(prob) 
     float prob; 
{ 
  float randomperc(); 
  
  if(randomperc() <= prob) 
    return(1); 
  else 
    return(0); 
} 

// initialization routine for randomnormaldeviate 
 
initrandomnormaldeviate() 
{ 
  rndcalcflag = 1; 
} 

 //normal noise with specified mean & std dev: mu & sigma

double noise(mu ,sigma) 
     double mu, sigma; 
{ 
  double randomnormaldeviate(); 
  
  return((randomnormaldeviate()*sigma) + mu); 
} 

 //random normal deviate after ACM algorithm 267 / Box-Muller Method 

double randomnormaldeviate() 
{ 
  double sqrt(), log(), sin(), cos(); 
  float randomperc();  
  double t, rndx1; 
  
  if(rndcalcflag) 
    { 
      rndx1 = sqrt(- 2.0*log((double) randomperc())); 
      t = 6.2831853072 * (double) randomperc(); 
      rndx2 = sin(t); 
      rndcalcflag = 0; 
      return(rndx1 * cos(t)); 
    } 
  else 
    { 
      rndcalcflag = 1; 
      return(rndx2); 
    } 
} 
 
float randomperc() 
{ 
  jrand++;
  
  if(jrand >= 55) 
    { 
      jrand = 1; 
      advance_random(); 
    } 
  return((float) oldrand[jrand]); 
} 

 //Pick a random integer between low and high 
  
int rnd(low, high) 
     int low,high; 
{ 
  int i; 
  float randomperc(); 
  
  if(low >= high) 
    i = low; 
  else 
    { 
      i = (randomperc() * (high - low + 1)) + low; 
      if(i > high) i = high; 
    } 
  return(i); 
} 

 //real random number between specified limits 

float rndreal(lo ,hi) 
     float lo, hi; 
{ 
  return((randomperc() * (hi - lo)) + lo); 
} 

 //Get random off and running 

warmup_random(random_seed) 
     float random_seed; 
{ 
  int j1, ii; 
  double new_random, prev_random; 
  
  oldrand[54] = random_seed; 
  new_random = 0.000000001; 
  prev_random = random_seed; 
  
  for(j1 = 1 ; j1 <= 54; j1++) 
    { 
      ii = (21*j1)%54; 
      oldrand[ii] = new_random; 
      new_random = prev_random-new_random; 
      if(new_random<0.0) new_random = new_random + 1.0; 
      prev_random = oldrand[ii]; 
    } 
  
  advance_random(); 
  advance_random(); 
  advance_random(); 
  
  jrand = 0; 
} 


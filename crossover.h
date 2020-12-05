/*This is the file for formulating the crossover process*/
void crossover(population *new_pop_ptr,population *mate_pop_ptr) ;

void crossover(population *new_pop_ptr,population *mate_pop_ptr)
{
  int i,j,k,n,p1,p2;
  int rnd1[optima],rnd2[optima],flag,tmp[2*optima];
  float rnd_no;

  new_pop_ptr->ind_ptr = &(new_pop_ptr->ind[0]);
  mate_pop_ptr->ind_ptr = &(mate_pop_ptr->ind[0]);

  for(i=0, n=0; i < popsize/2; i++, n=n+2)
    {
      flag = 0;
      p1 = rnd(0, (popsize-1));
      p2 = rnd(0, (popsize-1));

      while(p1 == p2)
	p2 = rnd(0, (popsize-1));

      for (j=0; j < chrom; j++)
	{
	  new_pop_ptr->ind[n].genes[j] = 0;
	  new_pop_ptr->ind[n+1].genes[j] = 0;

	  if((mate_pop_ptr->ind[p1].genes[j]==1) || (mate_pop_ptr->ind[p2].genes[j]==1))
	    {
	      tmp[flag] = j;
	      flag++;
	    }
	}

	for(j=0; j < flag; j++) {
	  	rnd_no =  randomperc();
		if(rnd_no < 0.5){
	  		new_pop_ptr->ind[n].genes[tmp[j]] = mate_pop_ptr->ind[p1].genes[tmp[j]];
			new_pop_ptr->ind[n+1].genes[tmp[j]] = mate_pop_ptr->ind[p2].genes[tmp[j]];
			}
		  else{
		    new_pop_ptr->ind[n].genes[tmp[j]] = mate_pop_ptr->ind[p2].genes[tmp[j]];
	  		new_pop_ptr->ind[n+1].genes[tmp[j]] = mate_pop_ptr->ind[p1].genes[tmp[j]];
			}
	}
 }
}





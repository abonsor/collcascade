
/* Collisional evolution */

#include <stdio.h>
#include <math.h>
#include<stdlib.h>

double sum_array(double a[], int num_elements)
{
  int i;
  double sum=0;
   for (i=1; i<num_elements; i++)
   {
	 sum = sum + a[i];
   }
   return(sum);
}

double round_to_digits(double value, int digits)
{
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;   
}


char* itoa(int i, char b[]){
    char const digit[] = "0123456789";
    char* p = b;
    if(i<0){
        *p++ = '-';
        i *= -1;
    }
    int shifter = i;
    do{ //Move to where representation ends
        ++p;
        shifter = shifter/10;
    }while(shifter);
    *p = '\0';
    do{ //Move back, inserting digits as u go
        *--p = digit[i%10];
        i = i/10;
    }while(i);
    return b;
}


int main(void)
{
    
    /*Variables*/
    double delta_t;
    double r_in;
    double r_out;
    double imax;
    double ecc;
    double au=1.5e11;
    double vrel;
    double delta;

    
    double alpha_p;
    double alpha;
    double mtot_0;
    double dbl, dmax;
    
    double mearth=5.97360e+24;
    double grav=6.67259e-11;
    double msun=1.98900e+30;

    double ntime1;
    int ntime;
    int itime, i, imk, jmi, nmax, ifv, id;
    int outputinterval;
    int flag_zero;
    char * word;
    int nbin;
    double fv_0;
    //start from the beginning or read file?
    int start;  // solids only? 
    int time_s=0, time_v=0, time_gas=0;
    int nf=1e2, diam_scrap, index_1000km;

    float x;

    
    double y2, y1, r, y1a, y1b;
    double fvarray[nf], n_k;
    double rho_s=3e3;// density of solids kg/m3
    double	rho_v=1000; // density of volatiles kg/m3
    double grain_a=0.01 ; //grain size mm // defunct parameter

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Open File
    FILE *ifp, *ofp,  *ins,  *ofparam, *ofins,  *ofer, *ofchi, *forigin_read, *forigin_write,  *ofdump, *ofratec;
    
    char *mode = "r";
    char outputFilename[] = "collouts.dat";
    char originfile_read[]="origin1000xxxxxx.dat",originfile_write[]="origin1000xxxxxx.dat";
    char infilenames[] = "in_s.dat";
    char inparam[]="inparam.in", erfile[]="errors.dat";
    char initials[]="initials.dat";


/* Open output files*/
    ofp = fopen(outputFilename, "w");
    ofparam=fopen(inparam, "r");
    ofer=fopen(erfile, "w");
    ofratec = fopen("ratec.dat", "w");
     
    if (ofp == NULL ||ofparam == NULL || ofratec == NULL || ofer == NULL)
            {   
              printf("Error! Could not open file\n"); 
              exit(-1); // must include stdlib.h 
            }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    //read in the parameters from inparam.in Note that the format must be exactly as specified here !//
    
    fscanf(ofparam, "%s\t %d\t %s\n",&word, &start, &word);
    if (start == 1){ printf(" Starting from the beginining\n");} else {printf(" starting from input files\n");}
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &delta_t,&word);
    printf("dt %e \n", delta_t);
    //rin
    fscanf(ofparam, "%s\t %lf\t %s\n ",&word, &r_in,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &r_out,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &imax,&word);
    printf("r_in %e \t rout %e \t imax %e \n", r_in, r_out, imax);
    fscanf(ofparam, "%s\t %lf\n",&word, &ecc);
    printf("ecc %e \n", ecc);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &vrel,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &delta,&word);
    printf("delta %e\t  \n", delta);
    fscanf(ofparam, "%s\t %lf\n",&word, &alpha_p);
    printf("  %e\t  \n",  alpha_p);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &mtot_0,&word);

    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &rho_s,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &dmax,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &dbl,&word);
    fscanf(ofparam, "%s\t %lf\n",&word, &ntime1);
    printf("Ntime %e\n", ntime1);
    fscanf(ofparam, "%s\t %d\n",&word, &outputinterval);
    printf("Outputinterval %d\n", outputinterval);
    printf(" ntime %e, outputinterval %d\n", ntime1, outputinterval);

    fscanf(ofparam, "%s\t %lf %s %s\n",&word, &grain_a, &word, &word);
    printf(" %e\n", grain_a);
    printf(" ntime %e, outputinterval %d\n", ntime1, outputinterval);

    fscanf(ofparam, "%s\t %lf",&word, &fv_0);
    printf(" ntime %e, outputinterval %d\n", ntime1, outputinterval);
    printf(" ntime %e, outputinterval %d\n", ntime1, outputinterval);

    fclose(ofparam);

//***************************************************/
 
    printf("Delta t %e, Rin %e, Rout %e\n", delta_t, r_in, r_out);
    printf("Imax %e, Ecc %e, vrel %e\n", imax, ecc,vrel);
    printf("delta %e, alpha_p %e, mtot0 %e\n", delta, alpha_p, mtot_0);
    printf("rho  s %e, v %e\n",  rho_s, rho_v);
    printf("dmax %e\t dbl %e\t ntime %e\t outputinterval %d\n", dmax, dbl, ntime1, outputinterval);

//***************************************************/


 /* Find the size of the diameter array nbin*/
    i=0;
    x=dmax;
    while (x>dbl)
        {
            x=pow(1-delta, 1./3.)*x;
            i=i+1;
        }

    nbin=i;
    printf(" \n Nbin: %d\t,  %e\t DBL %e\n", nbin, x, dbl);
/* determine norigin at 10m (the code does not trace origin in bodies less than 10m*/
    int norigin;
    i=0;
    x=dmax;
    while (x>10.0)
      {
        x=pow(1-delta, 1./3.)*x;
        i=i+1;
      }

    norigin=i;
    printf(" \n Norigin: %d\n", norigin);
    
    //***************************************************/
      /// declare variable arrays which are a function of nbin

    double const1,q, rho_k, fv[nbin],rmid,mass_i[nbin], mass_s_k[nbin], diam_i[nbin], mtot_k[nbin] ;
    double rate_c[nbin],   mtot_initial[nbin],mtot_initial1[nbin];
    double pik, vol, mdots_k,   m_all, sum_mass;
    printf(" \n Nbin: %d\t,  %e\t DBL %e\n", nbin, x, dbl);
    double sum_origin, sum_origin_bin, mtotk_new;
 
    
//// These two arrays trace where the mass in each bin 'originated' at each timestep. Origin_old is used to calculate the origin at the next timestep origin_new and updated at each timestep.
    double *origin_new = malloc(sizeof(*origin_new) * nbin * nbin);
    if (!origin_new)
        {
            printf("couldn't malloc origin_new\n");
            exit(1);
        }
    
    double *origin_old = malloc(sizeof(*origin_old) * nbin * nbin);
    if (!origin_old)
        {
            printf("couldn't malloc origin_old\n");
            exit(1);
        }
  
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// track where the material originated first index is the current bin and second index is the original bin! origin[k][]/mtot_k[k] is the fractional mass in the kth bin from each other bin! 

   printf(" origin declared\n Nbin: %d\t,  %e\t DBL %e\n", nbin, x, dbl);
    
    /* assign initial array */
    diam_i[1]=dmax;
    mass_i[1]=rho_s*pow(diam_i[1],3.)*M_PI/6.;
    mass_s_k[1]=mass_i[1];

    printf("Mk : %e\t %e\t %e\t %e \n", mass_i[1], mass_s_k[1], diam_i[1], fv_0);
		      
    for (i=2;i<nbin;i++)
		{
            mass_i[i]=mass_i[i-1]*(1-delta);
            mass_s_k[i]=mass_i[i];
            diam_i[i]=pow(mass_i[i]*6./(M_PI*rho_s),1./3.);
            if (diam_i[i] > 1000. ) {index_1000km=i;}
            nmax=i;
		}
      
    printf("Index 1000km : %d\n", index_1000km);
    printf("Mk : %e\t %e\t %e\t %e \n", mass_i[1], mass_s_k[1], diam_i[1], fv_0);

    //***************************************************/
      /// mass_s_i  (mass in solids) is fixed, whilst mass_i (total mass) varies as volatiles are lost
  
    /*Primordial size distribution with a constant power law*/
           
    alpha=(alpha_p+2)/3; /// changing 26th April 2022 from alpha = alpha_p-1/3;
    fprintf(ofp, "%e \t %e \t ",delta_t, alpha);
    for (i=1;i<nbin;i++)
        {
          mtot_initial[i]=pow(mass_i[i],(2-alpha));
        }
    const1=mtot_0*mearth/sum_array(mtot_initial,nbin);

    for (i=1;i<nbin;i++)
        {
            mtot_initial1[i]=mtot_initial[i]*const1;
        }
    m_all=sum_array(mtot_initial1, nbin);
    printf("M all: %e  %e\n", m_all, mtot_0);
    printf("Mtot last bin%e\n",mtot_initial1[nbin-1]);

      //// check that there is not more mass in the first bin than the mass of one object! 
    if (mtot_initial1[1] < mass_i[1])
        {
            printf("ERROR %e %e\n", mtot_initial1[1], mass_i[1]);
        }

    //***************************************************/

      /* QD star and QS star Dispersal and shattering thresholds */
      
      float Qa = 620 ;/*$J kg$^{−1}$, */
      float a = 0.3;
      float Qb = 5.6e-3;/*$J kg$^{−1}$ and $*/
      float b = 1.5 ;
      double qdstar[nbin], qsstar[nbin], mck[nbin], mlr[nbin];
      int irk[nbin], ick[nbin], j, k, p;

      for (i=1;i<nbin;i++) 
        {
            qdstar[i]=Qa *pow(diam_i[i],-a) + Qb *pow(diam_i[i],b);
            qsstar[i]=Qa *pow(diam_i[i],-a);
            mck[i] = ( 2*qdstar[i]/pow(vrel,2) )*mass_i[i];
            for (j=1; mass_i[j]>mck[i]; j++)
                {
                    ick[i]=j;
                }
            if (mass_i[1] <mck[i])
                {
                    ick[i]=0;
                }
            mlr[i] = ( 2*qsstar[i]/pow(vrel,2) )*mass_i[i];
            for (j=1; mass_i[j]>mlr[i]; j++)
                {
                    irk[i]=j;
                }
       }
     
    //***************************************************/
      /* Intrinsic collision probability*/
      rmid=r_in+(r_out-r_in)/2.;
      vrel=ecc*pow((grav*msun)/(rmid*au), 0.5);
      printf("ecc: %e\n", ecc);
      printf("vrel: %e\n", vrel);
      vol=8*M_PI*pow(rmid*au,3)*ecc*sin(imax)*(1+pow(ecc,2)/3.);
   
      pik=M_PI*vrel/vol;
      
      printf("IIII Pik: %e\n", pik);
      fprintf(ofp, "%e \t %e \t ",pik, rmid);

      /* Redistribution Function */
      float kminusi[nbin], f_s[nbin];
      for (i=1; i<nbin; i++)
        {
            kminusi[i]=i;
        }
  
    for (i=1; i<nbin; i++)
        {
            f_s[i]= pow((1-delta),(kminusi[i]*(2-alpha)))*delta*(2-alpha)* pow((1./2.),((alpha-2)));
        }
  

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (start==1) /* Starting from the beginning */
        {
            printf("Made it to origin");
   		      
            for (k=1; k<norigin; k++)

                {
                    sprintf(originfile_write, "1_origin%04d.dat", k);
                    forigin_write=fopen(originfile_write, "w");
                    sum_origin=0;
	
/* Initiate the array which defines where the mass originated origin[] - this array instead of being two dimensional origin[k][i] being the fraction of mass in the kth bin which originated from the ith bin, where D_i>D_k, instead we take origin[k*norigin +i] */
            /* Determine the size norigin, as there is no need to trace the origin of bodies smaller than D[norigin]*/
                    for (i = 1; i < norigin; i++)
                        {
			  
                            if (k==i)
                                {
                                    size_t offset = k * norigin + i;
                                    origin_new[offset]=1.0;
                                }
                            else
                                {
                                    size_t offset = k * norigin + i;
                                    origin_new[offset]=0;
                                }

                            fprintf(forigin_write,"%d \t %lf\n", i, origin_new[k*norigin+i]);
                            sum_origin=sum_origin+origin_new[k*norigin+i];
                        }

                    fclose(forigin_write);
                }

        }
        else  /* Starting from the output files start !=1*/
        {
            printf("Reading origin from latest file!");
            // all detail in the next loop
        }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //either use initial size distribution, or read from a file !
      if (start==1) 
            {
	  //initialise solids
                for (i=1; i<nbin; i++)
                    {
                        mtot_k[i]=mtot_initial1[i];
                    }

	  // open a file to write the initial conditions too...
                ofins = fopen(initials, "w");
        
                if (ofins == NULL )
                    {
                        printf("Error! Could not open file\n");
                        exit(-1); // must include stdlib.h
                    }

                for (i=1; i<nbin; i++)
                    {
                        fprintf(ofins, "%e \t ",mtot_k[i]);
                    }

                printf("Made it to origin");
                fclose(ofins);
            }
        else /* Read in files - starting from nth timestep*/
            {
                if (infilenames == NULL )
                    {
                        printf("Error! Could not open file\n");
                        exit(-1); // must include stdlib.h
                    }
	  /* Commands for making input files from old output to restart
    grep -A 1 "Timestep: 73320" collouts.dat > in_s.dat
	  grep -A 1 "Timestep: 73320" colloutv.dat > in_v.dat 
	  grep -A 1 Gas: collouts.dat > gas.dat
	  grep -B 1 "Timestep: 125000" gas.dat > in_gas.dat
cp collouts.dat collouts_old.dat
cp colloutv.dat colloutv_old.dat
Make sure you copy to old! 
	  */

	
                ins = fopen(infilenames, "r");
                fscanf(ins, "%s\t %d\n",&word, &time_s);
                printf("%c\t %d \n",word, time_s );
                printf("Starting at S: %d \t \n", time_s);
                fprintf(ofp, "\n %d\n   ",time_s);
                fclose(ofp);
                ofp = fopen(outputFilename, "w");

                for (i=1; i<nbin; i++)
                    {
                        fscanf(ins, "%lf\t", &mtot_k[i]);
                        fprintf(ofp, "%e \t ",mtot_k[i]);
                    }
	   
	    /////////////////////////////////////////////////////////////////////
	  //// need to read in correctly the origin from the latest file! 

                for (k=1; k<norigin; k++)
                    {
                        sprintf(originfile_read, "Final_origin%04d.dat", k);
                        sprintf(originfile_write, "origin%04d%04d.dat", k, time_s);
                        forigin_read=fopen(originfile_read, "r");
                        if (forigin_read == NULL )
                            {
                                printf("Error! Could not open file\n");
                                exit(-1); // must include stdlib.h
                            }
               
                        fscanf(forigin_read, "%*[^\n]\n");

                        sum_origin=0;

		      /// read in the origin file! 
                        for (i = 1; i < norigin; i++)
                            {
                                fscanf(forigin_read,"%lf\n", &origin_new[k*norigin+i]);
                                sum_origin=sum_origin+origin_new[k*norigin+i];
                            }
                        printf("\n SUM ORIGIN %d, %e\n", k, sum_origin);
                        fclose(forigin_read);

                    } /* k loop*/
    
            } /* else restart */
     

/***********************************************************/

      printf("About to start timestep loop\n");
      /* sum over each timestep*/
      for (itime=time_s; itime<ntime1; itime++)
	    {
	      //printf(" Timestep: %d\n   ", itime);
	      // write to dump files alternatively, such that the code can be restarted at any timestep
	      if (itime % 2)
                {
                    ofdump=fopen("1_dump_collouts.dat", "w");
                }
	      else
                {
                    ofdump=fopen("2_dump_collouts.dat", "w");
                }

	      fprintf(ofdump, "\n  Timestep: %d\n   ", itime);
	      if (itime % outputinterval == 0) 
                    {
                        fprintf(ofp, "\n  Timestep: %d\n   ", itime);
                        fprintf(ofratec, "\n Timestep: %d\n   ",itime);
                        printf("\n Printing to file:  \n Timestep: %d\n   ", itime);
                    }
	    
	      for (j=1; j<nbin; j++)
		    {
		      rate_c[j]=0;
		    }

	      for (j=1; j<nbin; j++)
		    {
		      /////////////////////////////////////////////////////////////////////////////////////////////////////
		    	      
		      /// check for problems with diameters
			  if (diam_i[j] <0)
			    {
			      //diam_i[j]=0;
			      printf("ERROR in DIAM NEGATIVE %e \n", diam_i[j]);
			    }
			  if (diam_i[j] != diam_i[j])
			    {
			      printf("ERROR in DIAM  NAN %e \t %e \t%e \t %e \t%e \n", diam_i[j], mtot_k[j], rho_k, n_k);
			      //diam_i[j]=0;
			    }

			/*************************************************/
		    
			/// CALCULATE RATE_C 
		 

                for (i=1;i<ick[j]; i++)
                    {
                        rate_c[j]=rate_c[j]+(mtot_k[i]/mass_s_k[i]) *pow((diam_i[i]+diam_i[j]),2.)*pik/4.;
                    }
		      /// check for problems with rates
                if (rate_c[j] <0)
                    {
                        rate_c[j]=0;
                        printf("ERROR in rate_c NEGATIVE %e \n", rate_c[j]);
                    }
                if (rate_c[j] != rate_c[j])
                    {
                        printf("ERROR in rate_c NAN %e \t %e \t%e \t %e \t%e \n", rate_c[j], mtot_k[j],  rho_k, n_k);
                        rate_c[j]=0;
                    }
			  
		
		  /*************************************************/ //Error spotted 29Nov2021 (add j<norigin)
                if (j<norigin)
                    {
                        for (i=0; i<norigin; i++)
                        {
                            origin_old[j*norigin+i]=origin_new[j*norigin+i];
                        }
                    }

            } // j

	    
	   
	  /*  sum over all bins */
	      for (k=1; k<nbin; k++)
                {
                    rho_k=rho_s;
                    n_k=mtot_k[k]/mass_s_k[k];
		  // diameters can change as volatiles lost
		  // adjust all diameters before calculations! 
		  
		     /* calculate loss from kth bin*/
	
		  		  	      
		  // mass gain - based on loss from other bins! (j)  */
                    imk=k-log(2)/delta;

	 
                    mdots_k=0;
		  

                    if (imk>0 )
                        {
                            for (i=1; i<imk; i++)
                                {
                                    p=k-i;
                                    mdots_k=mdots_k+f_s[p]*mtot_k[i]*rate_c[i];
                                }
                        }

		  // origin_old represents the bins 0 to nbin for kth

		  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    sum_mass=0;
		   //// WHERE DOES THIS MASS COME FROM?

		   /// addition to k from i: origin[k][i] = + mdots_k_array[i]*delta_t
		   // total loss is  - rate_c[k]*mtot_k[k]*deltat [ loose mass equally as a fraction of mass currently originating from the ith bin] 
		   // Store as fractional origin
		  //fractional originating from bin i   origin[k][i]/sum(origin[k][:])

	  
                    mtotk_new=mtot_k[k]+(mdots_k-rate_c[k]*mtot_k[k])*delta_t;
		          
                    if (k<norigin)
                        {
                            sum_origin=0;
                            for (i = 1; i < norigin; i++)
                                {
                                    sum_origin=sum_origin+origin_old[k*norigin+i];
			                        if( (k<i)& (origin_old[k*norigin+i]>0)) /// test if there is material originating from bins where it can't physically => only bins with >2M_k
                                        {
                                            fprintf(ofer, "\n READ Origin greater than 0 in bin with mass <2M_k (unphysical) %d %d %d %d", k, i, itime, k);
                                        }
                                }
		      
		      ///// do the kth bin first!
                            origin_new[k*norigin+k]=(origin_old[k*norigin+k]*mtot_k[k]- origin_old[k*norigin+k]*mtot_k[k]*rate_c[k]*delta_t)/mtotk_new;

		      /// for each bin i

                            if (imk>0)
                                {
                                    for (i = 1; i < imk; i++)
                                        {
                                            sum_origin_bin=0;
                                /// summ the origins from all bins j //
                                
                                            for (j = 1; j < imk; j++)
                                                {
                                                    p=k-j;
                                                // collect from each j bin which bin originated from ith bin into k !
                                                    sum_origin_bin=sum_origin_bin + origin_old[j*norigin+i]*f_s[p]*rate_c[j]*mtot_k[j]*delta_t;
                                                }
                                            origin_new[k*norigin+i]=(origin_old[k*norigin+i]*mtot_k[k] + sum_origin_bin- origin_old[k*norigin+i]*mtot_k[k]*rate_c[k]*delta_t)/mtotk_new;
                                
                                
				  /// make a resolution limit
                                            if (origin_new[k*norigin+i] < 1e-15)
                                                {
                                                    origin_new[k*norigin+i]=0;
                                                }
                                            sum_mass=sum_mass + f_s[p]*rate_c[i]*mtot_k[i]*delta_t;

                                        }
                                }// imk
	      
                            if (itime % outputinterval == 0)
                                {
                                    sprintf(originfile_write, "origin%04d.dat", k);
                                    forigin_write=fopen(originfile_write, "a");
                                    fprintf(forigin_write, "\n#Timestep:%d\n   ", itime);
                                    for (i=1; i<norigin; i++)
                                        {
                                            fprintf(forigin_write, "%e\n ", origin_new[k*norigin+i]);
                                        }
                                    fclose(forigin_write);
            
                                }// outputime
                        } // if k<norigin

 	    
		  //// ADD/LOOSE MASS in SOLIDS to kth bin
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		    mtot_k[k]=mtot_k[k]+(mdots_k-rate_c[k]*mtot_k[k])*delta_t;
            if (mtot_k <0) // check for mass below zero - unphysical
                {
                    fprintf(ofer, "Error: Solid mass k: %d mtot_k: %e itime: %d\n", k, mtot_k[k], itime);
                    mtot_k[k]=0;
                }
		   if (mtot_k[k] != mtot_k[k])
                {
                    mtot_k[k]=0;
                } //handle Nan
		  
		   if (itime % outputinterval == 0) 
                {
                    fprintf(ofp, "%e \t ",mtot_k[k]);
                    fprintf(ofratec,"%e \t ",rate_c[k]);
                    if (k==nbin-4)
                       {
                           printf("%e \n ",mtot_k[k]);
                       }
                }
		   fprintf(ofdump, "%e \t ",mtot_k[k]);
	     
            } //k

        m_all=sum_array(mtot_k, nbin);

        fclose(ofdump);
	    
	    } //itime

    free(origin_new);
    free(origin_old);
    fclose(ofp);
    fclose(ofer);
    fclose(ofratec);

      return 0;
}

//
//  initiation.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/4/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef initiation_hpp
#define initiation_hpp

#include <iostream>
#include <memory>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "random_uniform.hpp"
using namespace blitz;

void initiation(int N0, int N0r, int N0K, double uniup_r, double uniup_K, double unilow_r, double unilow_K, double muhatr, double muhatK, double sigmahatr, double sigmahatK, Array<float,2> &initial_array0, double mix_ratio_initial)
{
    double initial_r_growth_rate[N0r];
    double initial_K_growth_rate[N0K];
    Array<float,2> radom_number(1,N0,FortranArray<2>());
    radom_number=random_uniform(N0);
    double rangr = uniup_r - unilow_r;
    for (int x=1; x<=N0r; x++)
    {
        double rand1 = (radom_number(1,x)*rangr)+unilow_r;
        initial_r_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand1, sigmahatr) + muhatr;
    }
    
    radom_number=random_uniform(N0);
    
    double rangK = uniup_K - unilow_K;
    for (int x=1; x<=N0K; x++)
    {
        double rand2 = (radom_number(1,x)*rangK)+unilow_K;
        initial_K_growth_rate[x-1]=gsl_cdf_gaussian_Pinv(rand2, sigmahatK) + muhatK;
    }
    int xxx=1;
    int yyy=1;
    int zzz=N0*mix_ratio_initial;
    for (int i=1;i<=N0;i=i+1)
    {
        if(i<=zzz)
        {
            initial_array0(i,1)=initial_r_growth_rate[xxx-1];
            initial_array0(i,2)=0;
            initial_array0(i,3)=0;
            initial_array0(i,4)=1;
            initial_array0(i,5)=initial_array0(i,1);
            initial_array0(i,6)=i;
            initial_array0(i,7)=1;
            initial_array0(i,8)=0;
            initial_array0(i,9)=0;
            xxx=xxx+1;
        }
        else
        {
            initial_array0(i,1)=initial_K_growth_rate[yyy-1];
            initial_array0(i,2)=0;
            initial_array0(i,3)=0;
            initial_array0(i,4)=0;
            initial_array0(i,5)=initial_array0(i,1);
            initial_array0(i,6)=i;
            initial_array0(i,7)=1;
            initial_array0(i,8)=0;
            initial_array0(i,9)=0;
            yyy=yyy+1;
        }
    }
}







#endif /* initiation_hpp */

//
//  division.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/5/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef division_hpp
#define division_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "expected_growth_rate.hpp"
using namespace blitz;
void division(int i, Array<float,2> &initial_array, Array<float,2> temp1, Array<float,2> initial_array_temp,int carrying_capacity_r, int carrying_capacity_K, float alpha, float beta, int Kcell_numbers, int rcell_numbers, double death_time_range,double deltah, double max_growth_rate_r, double max_growth_rate_K, int C1)
{
    Range all = Range::all();
    const gsl_rng_type *T5;
    gsl_rng *r5;
    gsl_rng_env_setup();
    T5 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r5 = gsl_rng_alloc(T5);
    double growth_rate_density;
    if (initial_array(i,4)==1)//r-cell
    {
        growth_rate_density=expected_growth_rate(i, initial_array, carrying_capacity_r, Kcell_numbers, C1, alpha);
    }
    else //K-cell
    {
        growth_rate_density=expected_growth_rate(i, initial_array, carrying_capacity_K, rcell_numbers, C1, beta);
    }
    if (growth_rate_density>=0)
    {
        double X1=initial_array(i,1)*(1-0.05);
        double X2=initial_array(i,1)*(1+0.05);
        initial_array(i,1)=(X2-X1)*gsl_rng_uniform(r5)+X1;
        temp1(1,1)=(X2-X1)*gsl_rng_uniform(r5)+X1;
        initial_array(i,5)=growth_rate_density;
        temp1(1,5)=initial_array(i,5);
        temp1(1,6)=initial_array(i,6);
        temp1(1,7)=initial_array(i,7);
        temp1(1,8)=initial_array(i,8);
        temp1(1,9)=initial_array(i,9);
    }
    else
    {
        double probability_to_death=1/death_time_range;
        initial_array(i,8)=gsl_ran_geometric(r5,probability_to_death);
        temp1(1,7)=initial_array(i,7);
        initial_array(i,5)=growth_rate_density;
    }
    initial_array(i,2)=0;
    temp1(1,1)=initial_array(i,1);
    temp1(1,2)=initial_array(i,2);
    temp1(1,4)=initial_array(i,4);
    
    if (initial_array(i,4)==1)
    {
        if (initial_array(i,1) > max_growth_rate_r)
        {
            initial_array(i,1) = max_growth_rate_r;
        }
    }
    else if(initial_array(i,4)==0)
    {
        if(initial_array(i,1) > max_growth_rate_K)
        {
            initial_array(i,1) = max_growth_rate_K;
        }
    }
    if (temp1(1,4)==1)
    {
        if(temp1(1,1) > max_growth_rate_r)
        {
            temp1(1,1) = max_growth_rate_r;
        }
    }
    else if (temp1(1,4)==0)
    {
        if(temp1(1,1) > max_growth_rate_K)
        {
            temp1(1,1) = max_growth_rate_K;
        }
    }
    int current_size=initial_array.rows();
    initial_array_temp.resize(current_size+1,9);
    initial_array_temp=0;
    for (int rows=1;rows<=current_size+1;rows++)
    {
        if(rows<=current_size)
        {
            initial_array_temp(rows,all)=initial_array(rows,all);
        }
        else
        {
            initial_array_temp(rows,all)=temp1(1,all);
        }
    }
    initial_array.resize(current_size+1,9);
    initial_array=0;
    initial_array(all,all)=initial_array_temp(all,all);
    gsl_rng_free(r5);
}

#endif /* division_hpp */

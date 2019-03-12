//
//  death_dudgement.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/5/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef death_judgement_hpp
#define death_judgement_hpp

#include <stdio.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "expected_growth_rate.hpp"
using namespace blitz;
void death_judgement(Array<float,2> &initial_array, Array<float,2> initial_array_temp,int carrying_capacity_r, int carrying_capacity_K, float alpha, float beta, int Kcell_numbers, int rcell_numbers, double death_time_range,double deltah, double h)
{
    Range all = Range::all();
    int C0= initial_array.rows();
    const gsl_rng_type *T4;
    gsl_rng *r4;
    gsl_rng_env_setup();
    T4 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r4 = gsl_rng_alloc(T4);
    double expected_dividing_time=0;
    for (int i=1; i<=C0; i++)
    {
        if (initial_array(i,5)>=0)
        {
            if (initial_array(i,4)==1)//r-cell
            {
                initial_array(i,5)=expected_growth_rate(i, initial_array, carrying_capacity_r, Kcell_numbers, C0, alpha);
            }
            else //K-cell
            {
                initial_array(i,5)=expected_growth_rate(i, initial_array, carrying_capacity_K, rcell_numbers, C0, beta);
            }
            if (initial_array(i,5)<=0)
            {
                double probability_to_death=1/death_time_range;
                initial_array(i,8)=gsl_ran_geometric(r4,probability_to_death);
                initial_array(i,9)=initial_array(i,9)+deltah;
            }
            else
            {
                double expected_division_time=24/initial_array(i,5);
                double diving_time_range=0.1*expected_division_time;
                double undividing_time=0.9*expected_division_time;
                double probability_of_division=1/diving_time_range;
                expected_dividing_time=undividing_time+gsl_ran_geometric(r4,probability_of_division);
                initial_array(i,3)=expected_dividing_time;
            }
        }
        else
        {
            if(h==0)
            {
                double probability_to_death=1/death_time_range;
                initial_array(i,8)=gsl_ran_geometric(r4,probability_to_death);
            }
            else
            {
                double D_time_1=1.5*(24/initial_array(i,1));
                double D_time_2=0.9*initial_array(i,8);
                double D_time = 0;
                if (D_time_1<=D_time_2)
                {
                    D_time = D_time_1;
                }
                else
                {
                    D_time = D_time_2;
                }
                if (initial_array(i,9)<=D_time)///re-calculate the density growth rate
                {
                    if (initial_array(i,4)==1)//r-cell
                    {
                        initial_array(i,5)=expected_growth_rate(i, initial_array, carrying_capacity_r, Kcell_numbers, C0, alpha);
                    }
                    else //K-cell
                    {
                        initial_array(i,5)=expected_growth_rate(i, initial_array, carrying_capacity_K, rcell_numbers, C0, beta);
                    }
                    if (initial_array(i,5)<=0)
                    {
                        double probability_to_death=1/death_time_range;
                        initial_array(i,8)=gsl_ran_geometric(r4,probability_to_death);
                        initial_array(i,9)=initial_array(i,9)+deltah;
                    }
                    else
                    {
                        double expected_division_time=24/initial_array(i,5);
                        double diving_time_range=0.1*expected_division_time;
                        double undividing_time=0.9*expected_division_time;
                        double probability_of_division=1/diving_time_range;
                        expected_dividing_time=undividing_time+gsl_ran_geometric(r4,probability_of_division);
                        initial_array(i,3)=expected_dividing_time;
                    }
                }
                else
                {
                    if (initial_array(i,8)<=initial_array(i,9))
                    {
                        initial_array(i,7)=0;
                    }
                    else
                    {
                        initial_array(i,9)=initial_array(i,9)+deltah;
                    }
                }
            }
        }
    }
    int current_size=initial_array.rows();
    int sum =0;
    for (int CN=1; CN<=current_size; CN++)
    {
        if (initial_array(CN,7)==1)
        {
            sum=sum+1;
        }
    }
    initial_array_temp.resize(sum,9);
    initial_array_temp=0;
    int site1=1;
    for (int site=1; site<= current_size; site++)
    {
        if (initial_array(site,7)==1)
        {
            initial_array_temp(site1,all)=initial_array(site,all);
            site1++;
        }
    }
    initial_array.resize(sum,9);
    initial_array=0;
    initial_array(all,all)=initial_array_temp(all,all);
    gsl_rng_free(r4);
}

#endif /* death_judgement_hpp */

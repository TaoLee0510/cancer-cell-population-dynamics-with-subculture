//
//  density_growth_rate_calculation.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/5/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef density_growth_rate_calculation_hpp
#define density_growth_rate_calculation_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "expected_growth_rate.hpp"

using namespace blitz;
void density_growth_rate_calculation(Array<float,2> &initial_array, int carrying_capacity_r, int carrying_capacity_K, float alpha, float beta, int Kcell_numbers0, int rcell_numbers0)
{
    int C0= initial_array.rows();
    const gsl_rng_type *T3;
    gsl_rng *r3;
    gsl_rng_env_setup();
    T3 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r3 = gsl_rng_alloc(T3);
    double expected_dividing_time=0;
    for (int i=1; i<=C0; i++)
    {
        if (initial_array(i,4)==1)//r-cell
        {
            if (initial_array(i,5)>=0)
            {
                initial_array(i,5)=expected_growth_rate(i, initial_array, carrying_capacity_r, Kcell_numbers0, C0, alpha);
            }
        }
        else //K-cell
        {
            if (initial_array(i,5)>=0)
            {
                initial_array(i,5)=expected_growth_rate(i, initial_array, carrying_capacity_K, rcell_numbers0, C0, beta);
            }
        }
        double expected_division_time=24/initial_array(i,5);
        double diving_time_range=0.1*expected_division_time;
        double undividing_time=0.9*expected_division_time;
        double probability_of_division=1/diving_time_range;
        expected_dividing_time=undividing_time+gsl_ran_geometric(r3, probability_of_division);
        if (initial_array(i,2) >= expected_dividing_time)
        {
            initial_array(i,2)=0;
            initial_array(i,3)=gsl_ran_geometric(r3, 0.5);
        }
        else
        {
            initial_array(i,3)=expected_dividing_time;
        }
    }
    gsl_rng_free(r3);
}




#endif /* density_growth_rate_calculation_hpp */

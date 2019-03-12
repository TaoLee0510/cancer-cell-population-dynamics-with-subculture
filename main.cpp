//
//  main.cpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/4/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

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
#include "initiation.hpp"
#include "density_growth_rate_calculation.hpp"
#include "sort_time_division.hpp"
#include "death_judgement.hpp"
#include "division.hpp"
#include "subculture.hpp"
#include "save_data.hpp"
using namespace std;
using namespace blitz;


/////////////////////////////////////////////////////Usage//////////////////////////////////////////////////////////
////>>>>>>>>>>>>>>>>>>>>>>>./Density\ depenent\ cell\ population\ dynamics\ simulation N0 R mix_ratio_initial alpha_begin alpha_max alpha_interval beta_begin beta_max beta_interval time_interval deltah passage passage_times;

int main(int argc, char *argv[])
{
//int main()
//{
    
    int N0=atoi(argv[1]);
    int R=atoi(argv[2]);
    double mix_ratio_initial=atof(argv[3]);
    float alpha_begin=atof(argv[4]);
    float alpha_max=atof(argv[5]);
    float alpha_interval=atof(argv[6]);
    float beta_begin=atof(argv[7]);
    float beta_max=atof(argv[8]);
    float beta_interval=atof(argv[9]);
    int time_interval=atof(argv[10]);
    double deltah=atof(argv[11]);
    int passage=atof(argv[12]);
    int subculture_time=atof(argv[13]);
    
  
//
//    int time_interval=72;
//    int passage=0;//////////1:yes; 0:no.
//    int N0=500;
//    int R=1000;
//    double mix_ratio_initial=0.9;
//    float alpha_begin=0;
//    float alpha_max=4;
//    float alpha_interval=0.5;
//    float beta_begin=0;
//    float beta_max=2;
//    float beta_interval=0.5;
//    double deltah=0.05;
    
//    int subculture_time=30;
    int carrying_capacity_r=N0*2;
    int carrying_capacity_K=1.147*carrying_capacity_r;
    double death_time_range=36;
   //////////////////////////////////////////////////////parameters//////////////////////////////////////////////////////////////////////////////////////
    double muhatr=1.1832;
    double sigmahatr=0.2441;
    double muhatK=0.6832;
    double sigmahatK=0.3764;
    double min_growth_rate_r=1.0722619;
    double min_growth_rate_K=0.33963482;
    double max_growth_rate_r=1.3171805;
    double max_growth_rate_K=0.99505180;
    ////////////////////////////////////////////////////Array definition/////////////////////////////////////////////////////////////////////
    Range all = Range::all();
    Array<float,2> initial_array0(N0,9,FortranArray<2>());
    initial_array0=0;
    Array<float,2> initial_array(N0,9,FortranArray<2>());
    initial_array=0;
    Array<float,2> initial_array_temp(1,9,FortranArray<2>());
    initial_array=0;
    Array<float,2> rk_ratio(R,subculture_time+1,FortranArray<2>());
    rk_ratio=0;
    Array<float,2> total_cell_number(1,subculture_time,FortranArray<2>());
    total_cell_number=0;
    Array<float,2> temp1(1,9,FortranArray<2>());
    temp1=0;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int N0r=N0*mix_ratio_initial;
    int N0K=N0-N0r;
    double unilow_r=gsl_cdf_gaussian_P(min_growth_rate_r-muhatr, sigmahatr );
    double uniup_r=gsl_cdf_gaussian_P(max_growth_rate_r-muhatr, sigmahatr );
    double unilow_K=gsl_cdf_gaussian_P(min_growth_rate_K-muhatK, sigmahatK );
    double uniup_K=gsl_cdf_gaussian_P(max_growth_rate_K-muhatK, sigmahatK );
    for (float beta=beta_begin; beta<=beta_max; beta=beta+beta_interval)
    {
        for (float alpha=alpha_begin; alpha<=alpha_max; alpha=alpha+alpha_interval)
        {
            double h=0;
            int H=0;
            if (passage==1)
            {
                for (int replicates=1;replicates<=R;replicates=replicates+1)
                {
                    initiation(N0, N0r, N0K, uniup_r, uniup_K, unilow_r, unilow_K, muhatr, muhatK, sigmahatr, sigmahatK, initial_array0, mix_ratio_initial);
                    initial_array(all,all)=initial_array0(all,all);//////////////////////////////////initiation//////////////////////////////////////////////
                    for (int SBT=1;SBT<=subculture_time;SBT++)
                    {
                        h=0;
                        int C0=initial_array.rows();
                        int rcell_numbers0=0;
                        for (int i=1;i<=C0;i++)
                        {
                            rcell_numbers0=rcell_numbers0+initial_array(i,4);
                        }
                        int Kcell_numbers0=C0-rcell_numbers0;
                        density_growth_rate_calculation(initial_array, carrying_capacity_r, carrying_capacity_K, alpha, beta, Kcell_numbers0, rcell_numbers0);//////////////////////////density growth rate calculation///////////////////
                        for (H=1;H<1000000;H=H+1)
                        {
                            int C=initial_array.rows();
                            int rcell_numbers=0;
                            for (int i=1;i<=C;i++)
                            {
                                rcell_numbers=rcell_numbers0+initial_array(i,4);
                            }
                            int Kcell_numbers=C-rcell_numbers;
                            if ( h>=time_interval+1)
                            {
                                break;
                            }
                            if (h!=0)
                            {
                                death_judgement(initial_array, initial_array_temp, carrying_capacity_r, carrying_capacity_K, alpha, beta, Kcell_numbers, rcell_numbers, death_time_range, deltah, h);
                            }
                            sort_time_division(initial_array,initial_array_temp);///////sort time division/////////
                            int C1=initial_array.rows();
                            for (int i=1; i<=C1; i++)
                            {
                                if (initial_array(i,5)>=0)///////alive
                                {
                                    if (initial_array(i,2) >= initial_array(i,3)) /////////division
                                    {
                                        division(i, initial_array, temp1, initial_array_temp,carrying_capacity_r, carrying_capacity_K, alpha, beta, Kcell_numbers, rcell_numbers, death_time_range, deltah, max_growth_rate_r, max_growth_rate_K, C1);
                                    }
                                    else ///////no division
                                    {
                                        initial_array(i,2)=initial_array(i,2)+deltah;
                                    }
                                }
                            }
                            h=h+deltah;
                        }
                        int C=initial_array.rows();
                        if (C>=N0)
                        {
                            subculture(C, initial_array, initial_array_temp,N0);////////////////////////subculture///////////////////
                        }
                        int C2=initial_array.rows();
                        int rcell_numbers=0;
                        for (int i=1;i<=C2;i++)
                        {
                            rcell_numbers=rcell_numbers+initial_array(i,4);
                        }
                        rk_ratio(replicates,SBT)=double(rcell_numbers)/double(C2);
                    }
                    cout << replicates << endl;
                    cout << rk_ratio(replicates,all) << endl;////////////////////////display///////////////////
                    initial_array.resize(N0,9);
                    initial_array=0;
                }
            }
            else if (passage==0)
            {
                initiation(N0, N0r, N0K, uniup_r, uniup_K, unilow_r, unilow_K, muhatr, muhatK, sigmahatr, sigmahatK, initial_array0, mix_ratio_initial);
                initial_array(all,all)=initial_array0(all,all);//////////////////////////////////initiation//////////////////////////////////////////////
                    h=0;
                    int C0=initial_array.rows();
                    int rcell_numbers0=0;
                    for (int i=1;i<=C0;i++)
                    {
                        rcell_numbers0=rcell_numbers0+initial_array(i,4);
                    }
                    int Kcell_numbers0=C0-rcell_numbers0;
                    density_growth_rate_calculation(initial_array, carrying_capacity_r, carrying_capacity_K, alpha, beta, Kcell_numbers0, rcell_numbers0);//////////////////////////density growth rate calculation///////////////////
                    for (int H=1;H<99999999;H=H+1)
                    {
                        int C=initial_array.rows();
                        int rcell_numbers=0;
                        for (int i=1;i<=C;i++)
                        {
                            rcell_numbers=rcell_numbers0+initial_array(i,4);
                        }
                        int Kcell_numbers=C-rcell_numbers;
                        if ( h>=time_interval)
                        {
                            break;
                        }
                        if (h!=0)
                        {
                            death_judgement(initial_array, initial_array_temp, carrying_capacity_r, carrying_capacity_K, alpha, beta, Kcell_numbers, rcell_numbers, death_time_range, deltah, h);
                        }
                        sort_time_division(initial_array,initial_array_temp);///////sort time division/////////
                        int C1=initial_array.rows();
                        for (int i=1; i<=C1; i++)
                        {
                            if (initial_array(i,5)>=0)///////alive
                            {
                                if (initial_array(i,2) >= initial_array(i,3)) /////////division
                                {
                                    division(i, initial_array, temp1, initial_array_temp,carrying_capacity_r, carrying_capacity_K, alpha, beta, Kcell_numbers, rcell_numbers, death_time_range, deltah, max_growth_rate_r, max_growth_rate_K, C1);
                                }
                                else ///////no division
                                {
                                    initial_array(i,2)=initial_array(i,2)+deltah;
                                }
                            }
                        }
                        save_data(alpha, beta, rk_ratio, R, subculture_time,passage,initial_array, H);//////////////////////////save data//////////////////
                        h=h+deltah;
                    }
            }
            if (passage==1)
            {
                save_data(alpha, beta, rk_ratio, R, subculture_time,passage,initial_array, H);
            }
        }
    }
    cout << "Dynamics prediction finished" << endl;
    return 0;
}

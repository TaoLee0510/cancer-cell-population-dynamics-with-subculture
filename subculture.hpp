//
//  subculture.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/5/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef subculture_hpp
#define subculture_hpp

#include <stdio.h>
#include <random>
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

using namespace blitz;
void subculture(int C, Array<float,2> &initial_array, Array<float,2> initial_array_temp,int N0)
{
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    Range all = Range::all();
    int index[C];
    for (int a=1;a<=C;a++)
    {
        index[a-1]=a;
    }
    shuffle(index,index+C,RNG);
    shuffle(index,index+C,RNG);
    shuffle(index,index+C,RNG);
    shuffle(index,index+C,RNG);
    initial_array_temp.resize(N0,9);
    for (int b=1;b<=N0;b++)
    {
        initial_array_temp(b,all)=initial_array(index[b-1],all);
    }
    initial_array.resize(N0,9);
    initial_array=0;
    initial_array(all,all)=initial_array_temp(all,all);
}

#endif /* subculture_hpp */

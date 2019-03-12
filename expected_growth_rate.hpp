//
//  expected_growth_rate.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/6/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef expected_growth_rate_hpp
#define expected_growth_rate_hpp

#include <stdio.h>



double expected_growth_rate(int i, Array<float,2> &initial_array, int carrying_capacity, int cell_numbers, int total_numbers, float competition_co)
{
    double growth_rate_inherent=initial_array(i,1);
    double V3;
    double D1=growth_rate_inherent/carrying_capacity;
    V3=D1*total_numbers;
    double V4=competition_co*cell_numbers*D1;
    double growth_rate_density=growth_rate_inherent-V3-V4;
    return growth_rate_density;
}



#endif /* expected_growth_rate_hpp */

//
//  sort_time_division.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/5/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef sort_time_division_hpp
#define sort_time_division_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace blitz;
struct element
{
    double data;
    size_t index;
};
int compare(const void *a, const void *b)
{
    return (*(const element*)a).data - (*(const element*)b).data;
}
void sort_time_division(Array<float, 2> &initial_array, Array<float, 2> initial_array_temp)
{
    Range all = Range::all();
    int C_16= initial_array.rows();
    double cell_array_16[C_16];
    for (int CN=0; CN<C_16; CN++)
    {
        cell_array_16[CN]=initial_array(CN+1,3);
    }
    struct element{
        double data;
        size_t index;
    }array[C_16];
    for (size_t ii = 0; ii < sizeof(cell_array_16)/sizeof(cell_array_16[0]); ii++) {
        array[ii].data = cell_array_16[ii]*10000;
        array[ii].index = ii + 1;
    }
    qsort(array, sizeof(cell_array_16)/sizeof(cell_array_16[0]), sizeof(array[0]), compare);
    int index_16[C_16];
    for (int in=0; in<C_16; in++)
    {
        index_16[in]=(int)array[in].index;
    }
    initial_array_temp.resize(C_16, 9);
    for (int CNx=0; CNx<C_16; CNx++)
    {
        int ind=index_16[CNx];
        initial_array_temp(CNx+1,all)=initial_array(ind,all);
    }
    initial_array(all,all)=0;
    initial_array(all,all)=initial_array_temp(all,all);
}

#endif /* sort_time_division_hpp */

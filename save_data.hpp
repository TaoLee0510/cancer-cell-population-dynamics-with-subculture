//
//  save_data.hpp
//  Density depenent cell population dynamics simulation
//
//  Created by Tao Lee on 5/5/18.
//  Copyright © 2018 Tao Lee. All rights reserved.
//

#ifndef save_data_hpp
#define save_data_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
void save_data(float alpha, float beta, Array<float,2> rk_ratio, int R, int subculture_time,int passage,Array<float,2> initial_array, int H)
{
    if (passage==1)
    {
        char filedir [100] = {'\0'};
        sprintf(filedir, "./rk_ratio_a_%.1f_b_%.1f.txt",alpha,beta);
        FILE * fid1;
        fid1=fopen (filedir,"w+");
        for (int rows=1;rows<=R;rows++)
        {
            for (int col=1;col<=subculture_time+1;col++)
            {
                if (col<subculture_time+1)
                {
                    fprintf(fid1,"%g\t",rk_ratio(rows,col));
                }
                else
                {
                    fprintf(fid1,"%g\n",rk_ratio(rows,col));
                }
            }
        }
        fclose(fid1);
    }
    else if(passage==0)
    {
        char dirname [100] = {'\0'};
        sprintf(dirname, "mkdir ./a_%.1f_b_%.1f",alpha,beta);
        system(dirname);
        int C=initial_array.rows();
        char filedir [100] = {'\0'};
        sprintf(filedir, "./a_%.1f_b_%.1f/cell_index_a_%.1f_b_%.1f_t_%d.txt",alpha,beta,alpha,beta,H);
        FILE * fid1;
        fid1=fopen (filedir,"w+");
        for (int rows=1;rows<=C;rows++)
        {
            fprintf(fid1,"%g\n",initial_array(rows,6));
        }
        fclose(fid1);
    }
}




#endif /* save_data_hpp */

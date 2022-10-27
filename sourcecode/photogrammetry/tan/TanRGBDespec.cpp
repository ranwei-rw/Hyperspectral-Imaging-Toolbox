#include "mex.h"

#include <stdio.h>
#include <string.h>
#include "zHighlightRemoval.h"
#include "zGlobal.h"

void mexFunction(int nlhs, 
                 mxArray *plhs[],
                 int nrhs,
                 const mxArray *prhs[])

{

    mwSize n_elements;
    mwSize number_of_dims;
    const mwSize *dim_array;
    double *img;
    
    int i_element;
    int element_per_band;
    int mrows, ncols, nbands;
    int j_col, i_row, k_band;
    int count;
    double *outputptr;
    

    zArray2D<s_rgbi> diff;
    zArray2D<s_rgbi> raw_img;

    /* Get the number of elements in the input argument */
    n_elements = mxGetNumberOfElements(prhs[0]);
    /* Get the number of dimensions in the input argument. Allocate the space for the return argument */
    number_of_dims = mxGetNumberOfDimensions(prhs[0]);
    /* Get the number of dimensions in the input argument. */
    dim_array = mxGetDimensions(prhs[0]);
    mrows  = dim_array[0];
    ncols  = dim_array[1];
    nbands = dim_array[2];
    mexPrintf("Specularity removal is started\r");
    /* Get the data */
    img = mxGetPr(prhs[0]);

    diff.zAllocate(mrows, ncols);

    element_per_band = mrows * ncols;

    raw_img.zAllocate(mrows, ncols);

    //	mexPrintf("Allocate");
    for(j_col = 0; j_col < ncols; j_col ++)
    {
        for(i_row = 0; i_row < mrows; i_row ++)
        {
            raw_img[i_row][j_col].r = float(img[i_row + j_col * mrows + 0 * element_per_band]);
            raw_img[i_row][j_col].g = float(img[i_row + j_col * mrows + 1 * element_per_band]);
            raw_img[i_row][j_col].b = float(img[i_row + j_col * mrows + 2 * element_per_band]);
        }
    }

    mexPrintf("Ready to be processd\r");
    zHighlightRemoval hr(raw_img, diff);
    mexPrintf("Specular has been removed\r");
    plhs[0]   = mxCreateNumericArray(number_of_dims, dim_array,  mxDOUBLE_CLASS, mxREAL);
    outputptr = mxGetPr(plhs[0]);
    mexPrintf("Array has been formed\r");

    for(j_col = 0; j_col < ncols; j_col ++)
    {
        for(i_row = 0; i_row < mrows; i_row ++)
        {
            outputptr[i_row + j_col * mrows + 0 * element_per_band] = double(diff[i_row][j_col].r);		
            outputptr[i_row + j_col * mrows + 1 * element_per_band] = double(diff[i_row][j_col].g);
            outputptr[i_row + j_col * mrows + 2 * element_per_band] = double(diff[i_row][j_col].b);
        }
    }

    mexPrintf("Specularity removal is finished\r");

}
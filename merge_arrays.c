#include "mex.h"
#include "matrix.h"
#include "math.h"

mxChar *get_suffix(double *SA, mwSize n, int index, const mxArray *cell)
{
    /*the great memory saver!*/
    mxChar *ref_str = mxGetChars(mxGetCell(cell, (mwSize) SA[index + n] - 1));
    /*truncating to suffix*/
    return ref_str + (unsigned short) SA[index] - 1;
}

mwSize suf_len(double *SA, mwSize n, int index, const mxArray *cell)
{
    mwSize L = mxGetNumberOfElements(mxGetCell(cell, (mwSize) SA[index + n] - 1));
    return L - (mwSize) SA[index];
}

void free_suffix(const char *suf, double *SA, mwSize n, int index)
{
    /*unused at the moment, since no mxArrayToString*/
    suf -= (unsigned short) SA[index] - 1;
    mxFree(suf);
}

int lexcmp(mxChar *s1, mwSize len1, mxChar *s2, mwSize len2)
{
   /*got this online.*/
   mwSize minlen;

   minlen = (len1 <= len2) ? len1 : len2;
   minlen++;  /*needed to get the while condition right, when using mxChar*/

   /*
    * Compare as many bytes as are in the smaller string.  If an
    *  inequality is found, return the difference of the differing
    *  bytes.
    */
   while (minlen--)
      if (*s1++ != *s2++)
         return ((*--s1 & 0377) - (*--s2 & 0377));

   /*
    * The strings compared equal for the length of the shorter.  Return
    *  the difference in their lengths.  (Thus, the strings must be of
    *  the same length to be equal.)
    */
   return (len1 - len2);
}

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    /*using merge-sort to combine the two given (sorted) suffix arrays.
      this can be improved by using binary_search() on the longer SA.
      convention is that the first SA argument is the larger/older one.
      this convention only affects the behaviour of the unique feature,
      where the references from the first SA are kept for duplicates.
      Alon Diament, September 2015.*/
    double *A1, *A2;
    const mxArray *ref;
    mwSize n1, n2, nref;
    double unique_SA;
    mwSize i1, i2, insert;
    mxChar *suf1, *suf2;
    mwSize len1, len2;
    double *outSA;
    mwSize nOut;
    mwSize r, c;
    int diff;
    double *tmp = NULL;
    
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("merge_array:nrhs",
                          "at least 3 arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("merge_array:nlhs",
                          "only a single output supported.");
    }
    
    if (!mxIsDouble(prhs[0]) ||
        mxIsComplex(prhs[0]) ||
        mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("merge_array:InvalidInputType",
                          "first input must be a Nx3 double matrix");
    }
    if (!mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[1]) ||
        mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("merge_array:InvalidInputType",
                          "second input must be a Nx3 double matrix");
    }
    if (!mxIsCell(prhs[2]) ||
        fmin(mxGetN(prhs[2]), mxGetM(prhs[2])) != 1) {
        mexErrMsgIdAndTxt("merge_array:InvalidInputType",
                          "third input must be a 1D cell array of strings");
    }
    
    if (nrhs < 4) {
        unique_SA = 1;
    } else {
        if ((mxGetM(prhs[3]) != 1) ||
            (mxGetN(prhs[3]) != 1)) {
            mexErrMsgIdAndTxt("merge_array:InvalidInputType",
                              "fourth (optional) input must be a scalar");
        }
        unique_SA = mxGetScalar(prhs[3]);
    }

    A1 = mxGetPr(prhs[0]);
    n1 = mxGetM(prhs[0]);
    A2 = mxGetPr(prhs[1]);
    n2 = mxGetM(prhs[1]);
    ref = prhs[2];
    nref = fmax(mxGetN(prhs[2]), mxGetM(prhs[2]));
    
    /*final test: verifying cell array content*/
    for (r=0; r<nref; r++) {
        if (!mxIsChar(mxGetCell(ref, r))) {
            mexErrMsgIdAndTxt("merge_array:InvalidInputType",
                              "third argument must be a cell array of string (error on index %d)", i1);
        }
    }
    
    /*GET TO WORK!
      create the output matrix */
    nOut = n1 + n2;
    plhs[0] = mxCreateDoubleMatrix(nOut, 3, mxREAL);
    /* get a pointer to the real data in the output matrix */
    outSA = mxGetPr(plhs[0]);
    
    i1 = 0;
    i2 = 0;
    insert = 0;

    suf1 = get_suffix(A1, n1, i1, ref);
    len1 = suf_len(A1, n1, i1, ref);
    suf2 = get_suffix(A2, n2, i2, ref);
    len2 = suf_len(A2, n2, i2, ref);
    
    while ((i1 < n1) && (i2 < n2))
    {
        diff = lexcmp(suf1, len1, suf2, len2);

    /*mexPrintf("%s ? %s = %d\n", suf1, suf2, diff);*/

        if (diff < 0) {
            for (c=0; c<3; c++)
                outSA[insert+c*nOut] = A1[i1+c*n1];
            i1++;
            insert++;
            suf1 = get_suffix(A1, n1, i1, ref);
            len1 = suf_len(A1, n1, i1, ref);
        } else if (diff > 0) {
            for (c=0; c<3; c++)
                outSA[insert+c*nOut] = A2[i2+c*n2];
            i2++;
            insert++;
            suf2 = get_suffix(A2, n2, i2, ref);
            len2 = suf_len(A2, n2, i2, ref);
        } else {
            if (unique_SA) {
                /*increment frequency count*/
                A1[i1+2*n1] += A2[i2+2*n2];
            } else {
                for (c=0; c<3; c++)
                    outSA[insert+c*nOut] = A2[i2+c*n2];
                insert++;
            }
            i2++;
            suf2 = get_suffix(A2, n2, i2, ref);
            len2 = suf_len(A2, n2, i2, ref);
        }
    }
    
    if (unique_SA) {
        /*possibly need to truncate outSA due to ignored duplicates*/
        if (insert < i1 + i2) {
            mwSize nT = nOut - (i1 + i2 - insert);
            plhs[0] = mxCreateDoubleMatrix(nT, 3, mxREAL);
            tmp = mxGetPr(plhs[0]);
            for (r=0; r<insert; r++){
                for (c=0; c<3; c++){
                    tmp[r+c*nT] = outSA[r+c*nOut];
                }
            }
            outSA = tmp;
            nOut = nT;
        }
    }
    
    /*one array is over, copy rest of suffixes*/
    for (r=i1; r<n1; r++) {
        for (c=0; c<3; c++)
            outSA[insert+c*nOut] = A1[r+c*n1];
        insert++;
    }
    for (r=i2; r<n2; r++) {
        for (c=0; c<3; c++)
            outSA[insert+c*nOut] = A2[r+c*n2];
        insert++;
    }
}

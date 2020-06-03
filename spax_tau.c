int spax_tau(argc, argv)
int argc; void *argv[];
{
    int k;
    int n2;
    float *sa = (float *) argv[0];
    int *ija = (int *) argv[1];
    float *x = (float *) argv[2];
    float *b = (float *) argv[3];
    int n = *(int *) argv[4]; /* n = EVEN x-dimension of array = y-dimension */
    int nm = *(int *) argv[5]; /* nm = total number of sa */

    n2 = n/2;
    for (k=0;k<nm;k++) b[(k%n)/n2+(k/n)/n2*2] += sa[k]*x[ija[k]];
    return(0);
}

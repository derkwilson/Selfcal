int spax_tau_4ff(argc, argv)
int argc; void *argv[];
{
    int k;
    float *sa = (float *) argv[0];
    int *ija = (int *) argv[1];
    float *x = (float *) argv[2];
    float *b = (float *) argv[3];
    int n = *(int *) argv[4]; /* n = P */
    int nm = *(int *) argv[5]; /* nm = total number of sa */

    for (k=0;k<nm;k++) b[(k%4)+(k/n)*4] += sa[k]*x[ija[k]];
    return(0);
}

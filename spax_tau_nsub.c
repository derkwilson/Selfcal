int spax_tau_nsub(argc, argv)
int argc; void *argv[];
{
    int k;
    float *sa = (float *) argv[0];
    int *ija = (int *) argv[1];
    float *x = (float *) argv[2];
    float *b = (float *) argv[3];
    int n = *(int *) argv[4]; /* n = P */
    int nm = *(int *) argv[5]; /* nm = total number of sa */
    int nsub = *(int *) argv[6]; /* nsub = number of subsets in Py */

    for (k=0;k<nm;k++) b[(k%n)/(n/nsub)+(k/n)*nsub] += sa[k]*x[ija[k]];
    return(0);
}

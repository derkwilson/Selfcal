int sptx_tau(argc, argv)
int argc; void *argv[];
{
    int k;
    int n2;
    float *sa = (float *) argv[0];
    int *ija = (int *) argv[1];
    float *x = (float *) argv[2];
    float *b = (float *) argv[3];
    int n = *(int *) argv[4];
    int nm = *(int *) argv[5];

    n2 = n/2;
    for (k=0;k<nm;k++) b[ija[k]] += sa[k]*x[(k%n)/n2+(k/n)/n2*2];
    return(0);
}

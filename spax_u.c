int spax_u(argc, argv)
int argc; void *argv[];
{
    int k;
    float *sa = (float *) argv[0];
    int *ija = (int *) argv[1];
    float *x = (float *) argv[2];
    float *b = (float *) argv[3];
    int n = *(int *) argv[4]; /* number of pixels */
    int nm = *(int *) argv[5];

    for (k=0;k<nm;k++) {
        b[k%n] += sa[k]*x[ija[k%n]+(k/n)*4];
        }
    return(0);
}

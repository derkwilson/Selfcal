int sptx_u_nsub(argc, argv)
int argc; void *argv[];
{
    int k;
    float *sa = (float *) argv[0];
    int *ija = (int *) argv[1];
    float *x = (float *) argv[2];
    float *b = (float *) argv[3];
    int n = *(int *) argv[4];
    int nm = *(int *) argv[5];
    int nsub = *(int *) argv[6];

    for (k=0;k<nm;k++) {
        b[ija[k%n]+(k/n)*nsub] += sa[k]*x[k%n];
        }
    return(0);
}
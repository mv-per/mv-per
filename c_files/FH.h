


typedef struct {
    double n_max, H_VS, A_1v;
    double P;
} pure_fh_vsm;




void FHPcalc(int n, point_t *point, const void *arg)
{
    const pure_fh_vsm *params = (const pure_fh_vsm *)arg;
    double theta = point->x[0]/params->n_max;
    double Pcalc = (params->n_max/params->H_VS*(theta/(1.0-theta)))*\
                    exp(params->A_1v*params->A_1v*theta/(1.0 + params->A_1v*theta));
    point->fx = fabs(params->P-Pcalc)/params->P;
}




double CalculateFH(double P, double param[])
{
    // Estimativa Inicial: Langmuir
    // double b = param[1]/param[0];
    double NCalcIni = CalculateLangmuir(&P, param[0], param[1])*.8;
    
    // Define a função de minimização
    // static auto MinFunc = [&](double x) {
    // return FHPcalc(fabs(x), P, param[0], param[1], param[2]);
    // };
    // double (*minfun_pointer)(double) = [](double x) { return MinFunc(x); };

    // Minimiza
    /* Brent */
    // double result = brent_zeroin(minfun_pointer, fabs(NCalcIni), 1e-16);
        // if (isinf(result) || isnan(result)){
    //     return 10000.;
    // }
    // else {return result;}
    
    point_t start_pure;
    start_pure.x = (double *) malloc(1 * sizeof(double));
    start_pure.x[0] = NCalcIni;

    pure_fh_vsm params_pure;
    params_pure.n_max = param[0];
    params_pure.H_VS = param[1];
    params_pure.A_1v = param[2];
    params_pure.P = P;

    // optimisation settings
    optimset_t optimset_pure;
    optimset_pure.tolx = 1e-16;    // tolerance on the simplex solutions coordinates
    optimset_pure.tolf = 1e-16;    // tolerance on the function value
    optimset_pure.max_iter = 1500; // maximum number of allowed iterations
    optimset_pure.max_eval = 1500; // maximum number of allowed function evaluations
    optimset_pure.verbose = 0;     // toggle verbose output during minimization
    point_t solution_pure;
    nelder_mead(1, &start_pure, &solution_pure, &FHPcalc, &params_pure, &optimset_pure);

    return solution_pure.x[0];
}

double MinimizeFH(double P, double n_exp, double param[])
{
    // Minimiza
    double n_calc = CalculateFH(P, param);
    return fabs(n_exp-n_calc);
}

double * CalculateLogGammaFH(double *xs, double **Aij, size_t* ncomp)
{
    /* Variable declaration */
    double sum_1;
    double *gamma = (double *) malloc (*ncomp+1 * sizeof(double));
    size_t i, j;

    /* calculate gamma */
    for (i = 0; i <= *ncomp; i++){
        sum_1 = 0.;
        for (j = 0; j <= *ncomp; j++){
            sum_1 += xs[j] / (Aij[i][j] + 1.0);
        }

        gamma[i] = exp((1.0 - 1.0/sum_1) - log(sum_1));
    }
    return gamma;
}

double ** GetAijMatrixFH(vector<vector<double>> param, size_t* ncomp)
{
    /* Variable Declaration */    
    double **Aij = (double **) malloc (*ncomp+1 * sizeof(double *));
    size_t i, j;
    for (i = 0; i<=*ncomp; i++){
        Aij[i] = (double *) malloc (*ncomp+1 * sizeof(double));
    }
    /* Obtain Aij matrix */

    for (i = 0; i < *ncomp; i++){
        for (j = 0; j < *ncomp; j++){
            if (i == j) {Aij[i][j] = 0.0;}
            else {
                Aij[i][j] = (param[i][2] + 1.0)/(param[j][2] + 1.0) - 1.0 ;}
        }
    }
    for (i = 0; i < *ncomp; i++){
        Aij[i][*ncomp] = param[i][2];
        Aij[*ncomp][i] = param[i][2];
    }
    Aij[*ncomp][*ncomp] = 0.0;
    return Aij;
}

void MinimizeFHVSMMixture(int n, point_t *point, const void *arg)
{
    const mix_vsm_params *params = (const mix_vsm_params *)arg;

    //Function computation
    /* Variable declaration */
    size_t i;
    size_t ncomp = params->y.size();
    double *x = (double *) malloc (ncomp * sizeof(double));
    double *x_s = (double *) malloc (ncomp+1 * sizeof(double));
    double piA_RT, AdsorbedFugacity, ObjectiveFunction, n_max_mix, theta, nm, logGammaVacancy;
    double **Aij;
    double *gammaF;

    /* Set adsorbed phase composition */
    for ( i = 0; i < ncomp; i++){
        x[i] = point->x[i];
        // printf("x%d = %f", i, x[i]);
    }
    nm = point->x[ncomp];
    // printf("nm = %f\n", nm);

    /* Calculate the n_max for the mixture */
    n_max_mix = 0.0;
    for (i = 0; i < ncomp; i++){
        n_max_mix += x[i]*params->param[i][0];
    }

    // Define Theta
    theta = nm / n_max_mix;

    for (i = 0; i < ncomp; i++){
        x_s[i] = x[i]*theta;
    }
    x_s[ncomp] = 1.0 - theta;

    Aij = GetAijMatrixFH(params->param, &ncomp);

    gammaF = CalculateLogGammaFH(x_s, Aij, &ncomp);
    logGammaVacancy = log(gammaF[ncomp]*x_s[ncomp]);
    // double *BulkFugacity = params->PropsFun(params->y, params->P);
    
    // printf("\n start");
    double *BulkFugacity = (double *) malloc (ncomp * sizeof(double));
    for (i = 0; i < ncomp; i++){
        BulkFugacity[i] = params->y[i]* params->P;
    }   
    
    ObjectiveFunction = 0.0;
    double sumx = 0.;
    // double CompositionObjectiveFunction = 0.0;
    for (i = 0; i < ncomp; i++){
        sumx += fabs(x[i]);
        piA_RT = ((params->param[i][0] - n_max_mix)/nm - 1.0) * logGammaVacancy;
        AdsorbedFugacity = gammaF[i]*x[i]*(nm/n_max_mix)*(params->param[i][0]/params->param[i][1]) * (exp(params->param[i][2])/(1.0 + params->param[i][2])) * exp(piA_RT);
        // printf("\n BFug = %f, AdFug = %f, fobj[%d] = %f\n", BulkFugacity[i], AdsorbedFugacity, i, fabs(BulkFugacity[i] - AdsorbedFugacity)/BulkFugacity[i]);
        // ObjectiveFunction += fabs(BulkFugacity[i] - AdsorbedFugacity)/BulkFugacity[i];
        ObjectiveFunction += fabs(BulkFugacity[i] - AdsorbedFugacity)/BulkFugacity[i] ;
        
    }
    // printf("fx = %f", ObjectiveFunction + fabs(sumx-1.0)/sumx);
    point->fx = ObjectiveFunction + 100./(double) n *fabs(sumx-1.0);
}




vector<double> CalculateMixtureFHVSM(double P, vector<double> y, vector<vector<double>> param)
{
    int n = y.size();

    double *ini = CalculateExtendedLangmuir(param, y, &P, &n);

    point_t start; // initial point
    start.x = (double *) malloc(n+1 * sizeof(double));
    for (int i = 0; i <= n; i++) {
        // printf("ini[d]: %f \t", ini[i]);
        start.x[i] = ini[i];
    }


    mix_vsm_params params;
    params.P = P;
    params.y = y;
    params.param = param;


    // ;
    // optimisation settings
    optimset_t optimset;
    optimset.tolx = 1e-16;    // tolerance on the simplex solutions coordinates
    optimset.tolf = 1e-16;    // tolerance on the function value
    optimset.max_iter = 1500; // maximum number of allowed iterations
    optimset.max_eval = 1500; // maximum number of allowed function evaluations
    optimset.verbose = 0;     // toggle verbose output during minimization

    // printf("%f, %d", optimset.tolf, optimset.max_iter);
    // printf("tamanho y: %d", params.y.size());
    point_t solution;
    nelder_mead(n+1, &start, &solution, &MinimizeFHVSMMixture, &params, &optimset);

    vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        result[i] = solution.x[i]*solution.x[n];
        printf("result[%d]: %f \t",i, solution.x[i]);
    }
    printf("result[%d]: %f \t",n, solution.x[n]);
    printf("f(x) = %e\n", solution.fx);

    return result;
}


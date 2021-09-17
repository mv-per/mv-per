


typedef struct {
    double n_max, H_VS, tao_iv, tao_vi, alpha_iv;
    double P;
} pure_nrtl_vsm;




void NRTLPcalc(int n, point_t *point, const void *arg)
{
    const pure_nrtl_vsm *params = (const pure_nrtl_vsm *)arg;
    double G_iv = exp(-params->alpha_iv * params->tao_iv);
    double G_vi = exp(-params->alpha_iv * params->tao_vi);
    double theta = point->x[0]/params->n_max;
    double Pcalc = (params->n_max/params->H_VS*(theta/(1.0-theta)))*\
                    exp(params->tao_iv*G_iv*G_iv/((G_iv-1.0)*pow((G_iv-1.0)*theta + 1.0, 2.0)) +\
                    params->tao_vi*G_vi*G_vi/((G_iv-1.0)*pow((G_iv-1.0)*theta - G_vi, 2.0)) -\
                    params->tao_vi*((G_iv-1.0) + params->tao_iv*G_iv*G_iv*(G_vi-1.0))/((G_iv-1.0)*(G_vi-1.0)));

    point->fx = fabs(params->P-Pcalc)/params->P;
}




double CalculateNRTL(double P, double param[])
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

    pure_nrtl_vsm params_pure;
    params_pure.n_max = param[0];
    params_pure.H_VS = param[1];
    params_pure.tao_iv = param[2];
    params_pure.tao_vi = param[3];
    params_pure.alpha_iv = param[4];
    params_pure.P = P;

    // optimisation settings
    optimset_t optimset_pure;
    optimset_pure.tolx = 1e-16;    // tolerance on the simplex solutions coordinates
    optimset_pure.tolf = 1e-16;    // tolerance on the function value
    optimset_pure.max_iter = 1500; // maximum number of allowed iterations
    optimset_pure.max_eval = 1500; // maximum number of allowed function evaluations
    optimset_pure.verbose = 0;     // toggle verbose output during minimization
    point_t solution_pure;
    nelder_mead(1, &start_pure, &solution_pure, &NRTLPcalc, &params_pure, &optimset_pure);

    return solution_pure.x[0];
}

double MinimizeNRTL(double P, double n_exp, double param[])
{
    // Minimiza
    double n_calc = CalculateNRTL(P, param);
    return fabs(n_exp-n_calc);
}

double * CalculateLogGammaNRTL(double *xs, double **taoij, double **Gij, size_t *ncomp)
{
    /* Variable declaration */
    double sum_1, sum_2, sum_3, sum_4, sum_5;
    double *gamma = (double *) malloc (*ncomp+1 * sizeof(double));
    size_t i, j, k;

    /* calculate gamma */
    for (i = 0; i <= *ncomp; i++){
        sum_1 = 0.;
        sum_2 = 0.;
        for (j = 0; j <= *ncomp; j++){
            sum_1 += taoij[j][i] * Gij[j][i] * xs[j];
            sum_2 += Gij[j][i] * xs[j];
        }

        sum_3 = 0.;
        for (j = 0; j <= *ncomp; j++){
            sum_4 = 0.;
            sum_5 = 0.;
            for (k = 0; k <= *ncomp; k++){
                sum_4 += Gij[k][j]*xs[k];
                sum_5 += xs[k]*taoij[k][j]*Gij[k][j];
            }
            sum_3 += xs[j]*Gij[i][j]/sum_4 * (taoij[i][j] - sum_5/sum_4);
        }

        gamma[i] = exp(sum_1/sum_2 + sum_3);
    }
    return gamma;
}

double ** GetGijMatrixNRTL(vector<vector<double>> param, size_t *ncomp, double ** taoij)
{
    /* Variable Declaration */    
    double **Gij = (double **) malloc (*ncomp+1 * sizeof(double *));
    size_t i, j;
    for (i = 0; i<=*ncomp; i++){
        Gij[i] = (double *) malloc (*ncomp+1 * sizeof(double));
        for (j = 0; j <= *ncomp; j++){
            if (i == j) { Gij[i][j] = 1.0;}
            else {Gij[i][j] = exp( - 0.3 * taoij[i][j]);}
        }
    }

    for (i = 0; i < *ncomp; i++){
        Gij[i][*ncomp] = exp(-param[i][4]*param[i][2]);
        Gij[*ncomp][i] = exp(-param[i][4]*param[i][3]);
    }
    return Gij;
}

double ** GetTaoijMatrixNRTL(vector<vector<double>> param, size_t *ncomp)
{
    /* Variable Declaration */    
    double **taoij = (double **) malloc (*ncomp+1 * sizeof(double *));

    size_t i, j;
    for (i = 0; i<=*ncomp; i++){
        taoij[i] = (double *) malloc (*ncomp+1 * sizeof(double));
            for (j = 0; j <= *ncomp; j++){
                if (i == j) {taoij[i][j] = 0.0;}
                else {taoij[i][j] = 1.0;}
            }
    }

    for (i = 0; i < *ncomp; i++){
        taoij[i][*ncomp] = param[i][2];
        taoij[*ncomp][i] = param[i][3];
    }

    return taoij;
}

void MinimizeNRTLVSMMixture(int n, point_t *point, const void *arg)
{
    const mix_vsm_params *params = (const mix_vsm_params *)arg;

    //Function computation
    /* Variable declaration */
    size_t i;
    size_t ncomp = params->y.size();
    double *x = (double *) malloc (ncomp * sizeof(double));
    double *x_s = (double *) malloc (ncomp+1 * sizeof(double));
    double piA_RT, AdsorbedFugacity, ObjectiveFunction, n_max_mix, theta, nm, logGammaVacancy;
    double *gammaF;

    /* Set adsorbed phase composition */
    for ( i = 0; i < ncomp; i++){
        x[i] = point->x[i];
    }
    nm = point->x[ncomp];

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

    double ** taoij = GetTaoijMatrixNRTL(params->param, &ncomp);
    double ** Gij = GetGijMatrixNRTL(params->param, &ncomp, taoij);

    gammaF = CalculateLogGammaNRTL(x_s, taoij, Gij, &ncomp);
    logGammaVacancy = log(gammaF[ncomp]*x_s[ncomp]);

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
        AdsorbedFugacity = gammaF[i]*x[i]*(nm/n_max_mix)*(params->param[i][0]/params->param[i][1]) * exp(-(params->param[i][3]+ Gij[i][ncomp]*params->param[i][2])) * exp(piA_RT);
        ObjectiveFunction += fabs(BulkFugacity[i] - AdsorbedFugacity)/BulkFugacity[i] ;
    }
    point->fx =  ObjectiveFunction + 10. * fabs(sumx-1.0);
}




vector<double> CalculateMixtureNRTLVSM(double P, vector<double> y, vector<vector<double>> param)
{
    int n = y.size();

    double *ini = CalculateExtendedLangmuir(param, y, &P, &n);
    ini[n] *= .8;
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
    nelder_mead(n+1, &start, &solution, &MinimizeNRTLVSMMixture, &params, &optimset);

    vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        result[i] = solution.x[i]*solution.x[n];
        printf("result[%d]: %f \t",i, solution.x[i]);
    }
    printf("result[%d]: %f \t",n, solution.x[n]);
    printf("f(x) = %e\n", solution.fx);

    return result;
}


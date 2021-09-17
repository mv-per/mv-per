

/*
	Include libraries
*/
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <malloc.h>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <float.h>
using namespace std;

/* 
	Include Header files
*/
#include "nelder_mead.h"
#include "GETBOUNDARIES.h"
#include "BRENT.h"
#include "LANGMUIR.h"
#include "FH.h"
#include "WILSON.h"
#include "NRTL.h"

typedef struct {
    vector<vector<double>> param;
    vector<double> y;
    double P;
} mix_vsm_params;

/*
	Defines the functions used for the VSM-Wilson
*/
struct WILSON {

	~WILSON() {}

	double PureNads(double P, vector<double> param)
	{
		/*
		Return the adsorbed quantity for a single values of P
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}

		double nads = CalculateWilson(P, param_new);
		return nads;
		
	}
	
	vector<double> PureNads2(vector<double> P, vector<double> param)
	{
		/*
		Return the adsorbed quantity for multiple values of P
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}
		vector<double> nads(P.size(), 0.0);
		#pragma omp parallel for
		for (size_t i = 0; i < P.size(); i++)
		{
        	nads[i] = CalculateWilson(P[i], param_new); 
		}
		return nads;
		
	}

	double PureNdiff(vector<double> P, vector<double> n_exp, vector<double> param)
	{
		/*
		Return the absolute difference between calculated and exoerimental data of multiple points
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}
		double diff = 0.0;
		#pragma omp parallel for
		for (std::size_t i = 0; i < P.size(); i++) {
			diff += MinimizeWilson(P[i], n_exp[i], param_new); 
			
		}
		return diff;
	}

	vector<double> MixtureNads(double P, vector<double> y, vector<vector<double>> params)
	{
		// size_t i;
		size_t ncomp = y.size();
		// printf("%d", ncomp);
		//vector<double> result = CalculateMixtureWilsonVSM(P, y, params);
		// printf("%f", result[0]);
		// vector<double> new_result(ncomp, 0.0);
		// printf("P = %.2f\n", P);
		// for (i = 0; i < ncomp; i++){
		// 	// printf("n[%d] = %.8f\n",i, result[i]*result[ncomp]);
		// 	new_result[i] = result[i]*result[ncomp];
		// }
		
		return CalculateMixtureWilsonVSM(P, y, params);
	}

	vector<vector<double>> MixtureNadsMultiP(vector<double> P, vector<vector<double>> y, vector<vector<double>> param)
	{
        /* Return the xi, nm */
		size_t i;
		size_t ncomp = param.size();
		size_t PressureSizeList = P.size();
		vector<vector<double>> result;

		#pragma omp parallel for
		for (i = 0; i < PressureSizeList; i++){
			result[i] = CalculateMixtureWilsonVSM(P[i], y[i], param);
			// for (j = 0; j < ncomp; j++){
			// 	new_result[i][j] = result[j]*result[ncomp];
			// }
		}
		
		return result;
	}
};



/*
	Defines the functions used for the VSM-Flory-Huggins
*/
struct FH {

	~FH() {}

	double PureNads(double P, vector<double> param)
	{
		/*
		Return the adsorbed quantity for a single values of P
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}

		double nads = CalculateFH(P, param_new);
		return nads;
		
	}
	
	vector<double> PureNads2(vector<double> P, vector<double> param)
	{
		/*
		Return the adsorbed quantity for multiple values of P
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}
		vector<double> nads(P.size(), 0.0);
		#pragma omp parallel for
		for (size_t i = 0; i < P.size(); i++)
		{
        	nads[i] = CalculateFH(P[i], param_new); 
		}
		return nads;
		
	}

	double PureNdiff(vector<double> P, vector<double> n_exp, vector<double> param)
	{
		/*
		Return the absolute difference between calculated and exoerimental data of multiple points
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}
		double diff = 0.0;
		#pragma omp parallel for
		for (std::size_t i = 0; i < P.size(); i++) {
			diff += MinimizeFH(P[i], n_exp[i], param_new); 
		}
		return diff;
	}

	vector<double> MixtureNads(double P, vector<double> y, vector<vector<double>> params)
	{
		// size_t i;
		size_t ncomp = y.size();
		// printf("%d", ncomp);
		vector<double> result = CalculateMixtureFHVSM(P, y, params);
		// printf("%f", result[0]);
		// vector<double> new_result(ncomp, 0.0);
		// printf("P = %.2f\n", P);
		// for (i = 0; i < ncomp; i++){
		// 	// printf("n[%d] = %.8f\n",i, result[i]*result[ncomp]);
		// 	new_result[i] = result[i]*result[ncomp];
		// }
		
		return result;
	}

	vector<vector<double>> MixtureNadsMultiP(vector<double> P, vector<vector<double>> y, vector<vector<double>> param)
	{
        /* Return the xi, nm */
		size_t i;
		size_t ncomp = param.size();
		size_t PressureSizeList = P.size();
		vector<vector<double>> result;

		// #pragma omp parallel for
		for (i = 0; i < PressureSizeList; i++){
			result[i] = CalculateMixtureFHVSM(P[i], y[i], param);
			// for (j = 0; j < ncomp; j++){
			// 	new_result[i][j] = result[j]*result[ncomp];
			// }
		}
		return result;
	}
};

struct NRTL {

	~NRTL() {}

	double PureNads(double P, vector<double> param)
	{
		/*
		Return the adsorbed quantity for a single values of P
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}

		double nads = CalculateNRTL(P, param_new);
		return nads;
		
	}
	
	vector<double> PureNads2(vector<double> P, vector<double> param)
	{
		/*
		Return the adsorbed quantity for multiple values of P
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}
		vector<double> nads(P.size(), 0.0);
		#pragma omp parallel for
		for (size_t i = 0; i < P.size(); i++)
		{
        	nads[i] = CalculateNRTL(P[i], param_new); 
		}
		return nads;
		
	}

	double PureNdiff(vector<double> P, vector<double> n_exp, vector<double> param)
	{
		/*
		Return the absolute difference between calculated and exoerimental data of multiple points
		*/
		double *param_new = (double *) malloc(param.size() * sizeof(double));
		for (size_t i = 0; i < param.size(); i++){
			param_new[i] = fabs(param[i]);
		}
		double diff = 0.0;
		#pragma omp parallel for
		for (std::size_t i = 0; i < P.size(); i++) {
			diff += MinimizeNRTL(P[i], n_exp[i], param_new); 
		}
		return diff;
	}

	vector<double> MixtureNads(double P, vector<double> y, vector<vector<double>> params)
	{
		// size_t i;
		size_t ncomp = y.size();
		// printf("%d", ncomp);
		vector<double> result = CalculateMixtureNRTLVSM(P, y, params);
		// printf("%f", result[0]);
		// vector<double> new_result(ncomp, 0.0);
		// printf("P = %.2f\n", P);
		// for (i = 0; i < ncomp; i++){
		// 	// printf("n[%d] = %.8f\n",i, result[i]*result[ncomp]);
		// 	new_result[i] = result[i]*result[ncomp];
		// }
		
		return result;
	}

	vector<vector<double>> MixtureNadsMultiP(vector<double> P, vector<vector<double>> y, vector<vector<double>> param)
	{
        /* Return the xi, nm */
		size_t i;
		size_t ncomp = param.size();
		size_t PressureSizeList = P.size();
		vector<vector<double>> result;

		// #pragma omp parallel for
		for (i = 0; i < PressureSizeList; i++){
			result[i] = CalculateMixtureNRTLVSM(P[i], y[i], param);
			// for (j = 0; j < ncomp; j++){
			// 	new_result[i][j] = result[j]*result[ncomp];
			// }
		}
		return result;
	}
};


namespace py = pybind11;

/*
	Defines the pybind functions
*/
PYBIND11_MODULE(VSM, m) {

	py::class_<WILSON>(m, "WILSON")
		.def(py::init<>())
		.def("PureNads", &WILSON::PureNads)
		.def("PureNads2", &WILSON::PureNads2)
		.def("PureNdiff", &WILSON::PureNdiff)
		.def("MixtureNads", &WILSON::MixtureNads)
		.def("MixtureNads2", &WILSON::MixtureNadsMultiP);

	py::class_<FH>(m, "FH")
		.def(py::init<>())
		.def("PureNads", &FH::PureNads)
		.def("PureNads2", &FH::PureNads2)
		.def("PureNdiff", &FH::PureNdiff)
		.def("MixtureNads", &FH::MixtureNads)
		.def("MixtureNads2", &FH::MixtureNadsMultiP);

	py::class_<NRTL>(m, "NRTL")
		.def(py::init<>())
		.def("PureNads", &NRTL::PureNads)
		.def("PureNads2", &NRTL::PureNads2)
		.def("PureNdiff", &NRTL::PureNdiff)
		.def("MixtureNads", &NRTL::MixtureNads)
		.def("MixtureNads2", &NRTL::MixtureNadsMultiP);


#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}


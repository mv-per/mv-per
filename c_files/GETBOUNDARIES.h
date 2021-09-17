

const double sqrt_2 = 1.414213562373095048801;
const double eps = 1.053671e-08;


// Find the sign change
double sign(double X1, double X2)
{
	if (X2 >= 0) return fabs(X1);
	else return -fabs(X1);
}

// Find Boundaries for the Brent Method - From MATLAB library
void getboundaries(double (*fun)(double), double &x, double* a, double* b){
    double fa, fb, fx, dx;
    //assign initial values
    fx = fun(x);
	*a = x; fa = fx; *b = x; fb = fx;

    	// find a and b
	if (x != 0.) dx = x/50.0;
	else dx = 1.0/50.0;

    size_t k = 0;
	// find change of sign
    get_a_and_b:
        while ((fa > 0.) == (fb > 0.)){
            dx = 1.414213562373095048801 * dx;
            *a = x - dx;
            fa = fun(*a);
            if ((fa>0.) != (fb>0.)) break;
            *b = x + dx;
            fb = fun(*b);
        }

    if ((fabs(fa) > .5 || fabs(fb) > .5) && k < 10){
        dx /= 2.;
        k++;
        goto get_a_and_b;
    }

}
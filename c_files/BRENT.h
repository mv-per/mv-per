

// Solves a single variable function 
double brent_zeroin(double (*fun)(double), double x, double tol)
{
	double a, b, c, d, e, fa, fb, fc, tol1, xm, p, q, r, s, fx;


	//Compute eps, the relative machine precision
	// eps = 1.1e-16;

	//If manually determination of eps
	// define_eps:
	// 	tol = tol/2.0;
	// 	tol1 = 1.0 + tol;
	// 	if (tol1 > 1.0) goto define_eps;
	// // std::cout << eps << std::endl;
	

	
	// Find the best initial estimates of negative and positive values
    getboundaries(fun, x, &a, &b);
	fx = fun(x);

	// Find the function values of the estimates
	fa = fun(a);
	fb = fun(b);

	//Begin step
	begin_20:
		c = a;
		fc = fa;
		d = b-a;
		e = d;
	begin_30:
		if (fabs(fc) >= fabs(fb)) goto convergence_test;
		a = b;
		b = c;
		c = a;
		fa = fb;
		fb = fc;
		fc = fa;

	convergence_test:
		tol1 = 2.0*eps*fabs(b) + 0.5*tol;
		xm = 0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) goto done;

	// Is bisection necessary
	if (fabs(e) < tol1) goto bisection;
	if (fabs(fa) <= fabs(fb)) goto bisection;

	// Is quadratic size_terpolation possible
	if (a != c) goto inverse_quadratic_size_terpolation;

	s = fb/fa;
	p = 2.0*xm*s;
	q = 1.0 - s;

	goto adjust_signs;

	inverse_quadratic_size_terpolation:
		q = fa/fc;
		r = fb/fc;
		s = fb/fa;
		p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
		q = (q - 1.0)*(r - 1.0)*(s - 1.0);

	adjust_signs:
		if (p > 0.0) {
			q = -q;
		}
		p = fabs(p);

	// Is size_terpolation acceptable
	if (2.0*p >= (3.0*xm*q - fabs(tol1*q)) || p >= fabs(0.5*e*q)) goto bisection;
	e = d;
	d = p/q;
	goto complete_step;

	bisection:
		d = xm;
		e = d;

	complete_step:
		a = b;
		fa = fb;
		if (fabs(d) > tol1) b = b + d;
		else b = b + sign(tol1, xm);
		fb = fun(b);
		if ((fb*(fc/fabs(fc))) > 0.) goto begin_20;
		goto begin_30;

	done:
		// printf("f(x) = %.8f \n", fabs(xm));
		double zeroin = b;
		if (zeroin == x) {
			tol /= 2.;
			goto begin_20;
		} 
		return zeroin;

}
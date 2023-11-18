__declspec(dllexport)
void dogleg(int n, double r[], int lr, double diag[], double qtb[],
	double delta, double x[], double wa1[], double wa2[]);
__declspec(dllexport)
double enorm(int n, double x[]);

__declspec(dllexport)
void fdjac1(void fcn(int n, double x[], double f[]),
	int n, double x[], double fvec[], double fjac[], int ldfjac,
	int ml, int mu, double epsfcn, double wa1[], double wa2[]);

__declspec(dllexport)
int fsolve(void fcn(int n, double x[], double fvec[]), int n,
	double x[], double fvec[], double tolerance);

__declspec(dllexport)
int hybrd(void fcn(int n, double x[], double fvec[]),
	int n, double x[], double fvec[], double xtol, int maxfev, int ml,
	int mu, double epsfcn, double diag[], int mode, double factor,
	int nfev, double fjac[], int ldfjac, double r[], int lr, double qtf[],
	double wa1[], double wa2[], double wa3[], double wa4[]);

__declspec(dllexport)
void qform(int m, int n, double q[], int ldq);

__declspec(dllexport)
void qrfac(int m, int n, double a[], int lda, bool pivot, int ipvt[],
	int lipvt, double rdiag[], double acnorm[]);

__declspec(dllexport)
void r1mpyq(int m, int n, double a[], int lda, double v[], double w[]);

__declspec(dllexport)
bool r1updt(int m, int n, double s[], int ls, double u[], double v[], double w[]);
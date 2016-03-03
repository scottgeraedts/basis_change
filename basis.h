#ifndef BASIS_H
#define BASIS_H
#include "utils.h"
#include <numeric>
#include <algorithm>

//#include<complex>

using namespace std;
extern"C"{
	void jacobi_theta_(int *n, complex<double> *z, complex<double> *tau, complex<double> *theta, int *sum);
	void weierstrass_sigma_( complex<double> *z, complex<double> *l1, complex<double> *l2, complex<double> *sigma,int *reduce);
}

class ModelTorus{
public:
	int NPhi,Ne,nStates,invNu;
	double Lx,Ly;
	complex<double> dsum;
	vector< complex<double> > ds;	
	Eigen::MatrixXd T1odd,T1even,Tinvert;
	complex<double> L1,L2;
	vector<int> states,family_size;
	vector<bool> invert_family;
	bool shrink;

	ModelTorus();
	void make_states();
	void basis_change();
	complex<double> landau_basis(int m);
	double my_basis(int m);
	complex<double> many_body_laughlin(vector<int> sites);
	complex<double> duncan_cfl(const vector<int> &sites);
	void make_translate_one();
	int find_translated(int, const vector<int> &heads);
	complex<double> duncan_test(int);
	void translation_eigenvectors();
	complex<double> single_determinant(int k, const vector<int> &eval_sites_element, const vector< complex<double> > &single_orbitals);
};
#endif

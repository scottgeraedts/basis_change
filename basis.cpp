#include "basis.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>

//reads an input file and sets some important variables
ModelTorus::ModelTorus(){
	ifstream infile;
	infile.open("params");
	infile>>NPhi;
	infile>>Ne;
	Lx=sqrt(2.*NPhi*M_PI);
	Ly=Lx;
	double tempreal,tempimag;
	for(int i=0;i<Ne;i++){
		infile>>tempreal>>tempimag;
		ds.push_back(complex<double>(tempreal*Lx/Ne,tempimag*Ly/Ne));
	}
	infile.close();
	L1=complex<double>(Lx,0); L2=complex<double>(0,Ly);
	shrink=false;
	make_states();	
	invNu=NPhi/Ne;
	dsum=accumulate(ds.begin(), ds.end(), complex<double>(0,0) );
	double dvar=0;
	for(int i=0;i<Ne;i++) dvar+=norm(dsum-ds[i]);
	cout<<"ds data: "<<dsum<<" "<<dvar/(2*M_PI*Ne)<<endl;
}
	
void ModelTorus::basis_change(){
	bool post_shrink=false;

	complex<double> Tphase=polar(1.,dsum.real()*Ly/(1.*invNu*NPhi));
	cout<<"T:"<<Tphase<<endl;

	Eigen::MatrixXcd changer=Eigen::MatrixXcd::Zero(nStates,nStates);
	Eigen::MatrixXcd to_determinant(Ne,Ne);
	Eigen::VectorXcd model_wf(nStates);

	vector< complex<double> > single_particle_orbitals(NPhi);
	for(int i=0;i<NPhi;i++) single_particle_orbitals[i]=landau_basis(i);
//	for(int i=-10;i<15;i++) cout<<i<<" "<<landau_basis(i)<<endl;

	//build matrix and vector
	vector<int> eval_sites_element, inverted_element;	
	int cycled_sites,old_sites;
	double sign;
	clock_t t=clock();
	//loop through all families
	for(int i=0;i<nStates;i++){
		if(shrink){
			cycled_sites=states[i];
			sign=1;
			model_wf(i)=0;
			eval_sites_element=bitset_to_pos(states[i],NPhi);	
			model_wf(i)=many_body_laughlin(eval_sites_element)*sqrt(family_size[i]);
		
			//loop through all members of this family
			for(int j=0;j<NPhi;j++){
				eval_sites_element=bitset_to_pos(cycled_sites,NPhi);	
				if(invert_family[i]) inverted_element=bitset_to_pos(invert_bits(cycled_sites,NPhi),NPhi); 
				//update matrix
				for(int k=0;k<nStates;k++){
					changer(i,k)+=sign*single_determinant(k,eval_sites_element,single_particle_orbitals)*sqrt((1.*family_size[k])/(1.*family_size[i]));
					if(invert_family[i]) changer(i,k)+=sign*single_determinant(k,inverted_element,single_particle_orbitals)*sqrt((1.*family_size[k])/(1.*family_size[i]));
				}
				//shift sites and figure out new sign
				sign*=-1;
				old_sites=cycled_sites;
				cycled_sites=cycle_bits(cycled_sites,NPhi);
				if(cycled_sites<old_sites) sign*=-1;
				if(cycled_sites==states[i]) break;
			}
		}else{//don;t shrink
//			cout<<(bitset<8>)states[i]<<endl;
			eval_sites_element=bitset_to_pos(states[i],NPhi);
//			model_wf(i)=duncan_cfl(eval_sites_element);
//			model_wf(i)=many_body_laughlin(eval_sites_element);
			for(int k=0;k<nStates;k++){
//				cout<<(bitset<5>)states[i]<<" "<<(bitset<5>)states[k]<<endl;
				changer(i,k)=single_determinant(k,eval_sites_element,single_particle_orbitals);
			}
		}			
	}

	for(int COM=0;COM<NPhi;COM++){
		cout<<"starting run with COM="<<COM<<endl;
		model_wf=Eigen::VectorXcd::Zero(nStates);
		for(int i=0;i<nStates;i++){
			eval_sites_element=bitset_to_pos(states[i],NPhi);
			if(accumulate(eval_sites_element.begin(),eval_sites_element.end(),0)%NPhi==COM) model_wf(i)=duncan_cfl(eval_sites_element);
			else model_wf(i)=0;
//			model_wf(i)=duncan_cfl(eval_sites_element);
		}
		
		if(model_wf.squaredNorm()>1e-10) model_wf/=sqrt(model_wf.squaredNorm());
		cout<<"model wf:"<<endl;
		for(int i=0;i<nStates;i++) cout<<model_wf(i)<<"\t\t"<<(bitset<8>)states[i]<<" "<<states[i]<<endl;
	//	cout<<endl<<changer.real()<<endl;
	
		cout<<"time to build matrix "<<(float)(clock()-t)/(CLOCKS_PER_SEC)<<endl;
	//	cout<<changer.real()<<endl;
	
		vector<int> head_indices,family_heads;
		int nHeads;
		if(post_shrink && !shrink){
			make_translate_one(); 
			cout<<"wave function translational symmetry: "<<(T1odd*model_wf+Tphase*model_wf).norm()<<endl;
			cout<<"matrix translation symmetry: "<<(T1odd*changer*T1odd.transpose()-lil_sign(Ne)*changer).norm()<<endl;
			//Tphase=complex<double>(1,0);
			//use translation symmetry to shrink A
			Eigen::EigenSolver<Eigen::MatrixXd> Tsols(T1odd);

			//if you haven't used symmetry yet to shrink the basis, do it now
			for(int i=0;i<nStates;i++){
				if( norm(Tsols.eigenvalues()(i)+Tphase)<1e-15  ){
					head_indices.push_back(i);
				}
			}
			nHeads=head_indices.size();
			Eigen::MatrixXcd shrinker=Eigen::MatrixXcd::Zero(nStates,head_indices.size());
			for(int i=0;i<(signed)head_indices.size();i++){
				shrinker.col(i)=Tsols.eigenvectors().col(head_indices[i]);			
				for(int j=0;j<nStates;j++){
					if(shrinker.col(i)(j).real()>1e-12){
						family_heads.push_back(j);
						break;
					}
				}
			}
	//		Eigen::MatrixXcd shrinker=translation_eigenvectors();
	//		cout<<Tsols.eigenvalues()<<endl;
			for(int i=0;i<(signed)family_heads.size();i++) cout<<(bitset<8>)states[family_heads[i]]<<endl;
	//		cout<<shrinker<<endl;
			changer=shrinker.conjugate().transpose()*changer*shrinker;
			model_wf=shrinker.conjugate().transpose()*model_wf;
		
			//inversion symmetry
	//		int newstate,newhead;
	//		vector<int> sites;
	//		Tinvert=Eigen::MatrixXd::Zero(nHeads,nHeads);
	//		for(int i=0;i<nHeads;i++){
	//			sites=bitset_to_pos(states[family_heads[i]],NPhi);
	//			newstate=0;
	//			for(int j=0;j<(signed)sites.size();j++) newstate = newstate | 1<<((NPhi-1-sites[j])%NPhi);
	//			newhead=find_translated(newstate,family_heads);
	//			Tinvert(i,abs(newhead))=1;
	//			if( newhead<0 && NPhi%2==0) Tinvert(i,abs(newhead))=-1;
	////			cout<<(bitset<8>)newstate<<" "<<(bitset<8>)states[family_heads[i]]<<" ";
	////			cout<<(bitset<8>)states[family_heads[abs(newhead)]]<<" "<<(newhead>0) - (newhead<0)<<endl;
	//		}
	//		//Tinvert(4,4)=1; Tinvert(7,7)=-1;
	//		cout<<"wave function inverssion symmetry: "<<(Tinvert*model_wf-model_wf).norm()<<endl;
	//		cout<<"matrix inversion symmetry: "<<(Tinvert*changer*Tinvert.transpose()-changer).norm()<<endl;
	////		cout<<changer.real()<<endl;
	////		cout<<model_wf<<endl;
	//		
	//		Eigen::EigenSolver<Eigen::MatrixXd> InvSol(Tinvert);
	//		head_indices.clear();
	//		//if you haven't used symmetry yet to shrink the basis, do it now
	//		for(int i=0;i<nHeads;i++){
	//			if( norm(InvSol.eigenvalues()(i)-complex<double>(1,0))<1e-15  ){
	//				head_indices.push_back(i);
	//			}
	//		}
	//		vector<int> old_family_heads=family_heads;
	//		family_heads.clear();
	//		Eigen::MatrixXcd shrinker2=Eigen::MatrixXcd::Zero(nHeads,head_indices.size());
	//		for(int i=0;i<(signed)head_indices.size();i++){
	//			shrinker2.col(i)=InvSol.eigenvectors().col(head_indices[i]);			
	//			for(int j=0;j<nHeads;j++){
	//				if(abs(shrinker2.col(i)(j))>1e-12){
	//					family_heads.push_back(old_family_heads[j]);
	//					break;
	//				}
	//			}
	//		}
	////		Eigen::MatrixXcd shrinker=translation_eigenvectors();
	////		cout<<Tsols.eigenvalues()<<endl;
	//		for(int i=0;i<(signed)family_heads.size();i++) cout<<(bitset<8>)states[family_heads[i]]<<endl;
	////		cout<<shrinker<<endl;
	//		changer=shrinker2.conjugate().transpose()*changer*shrinker2;
	//		model_wf=shrinker2.conjugate().transpose()*model_wf;
		}
	
		if(!post_shrink && !shrink){	//some symmetry tests
			make_translate_one();
			cout<<"wave function symmetry: "<<(T1odd*model_wf+Tphase*model_wf).norm()<<endl;
	//		cout<<"wave function symmetry: "<<(Tinvert*model_wf-model_wf).norm()<<endl;
			Eigen::MatrixXcd tmw=Tinvert*model_wf;
			for(int i=0;i<nStates;i++)
				if (abs(tmw(i)-model_wf(i))>1e-10) cout<<model_wf(i)<<" "<<tmw(i)<<" "<<(bitset<8>)states[i]<<endl;
			cout<<"matrix symmetry: "<<(T1odd*changer*T1odd.transpose()-lil_sign(Ne)*changer).norm()<<endl;
	//		cout<<"matrix symmetry: "<<(Tinvert*changer*Tinvert.transpose()-changer).norm()<<endl;
		}
	
	//	//solve changer*x=model_wf	
		t=clock();
		Eigen::FullPivLU<Eigen::MatrixXcd> lu(changer);
		cout<<"determinant"<<lu.determinant()<<endl;
	//	cout<<"changer is invertible? "<<lu.isInvertible()<<endl;
	//	Eigen::VectorXcd out=lu.solve(model_wf);
		Eigen::VectorXcd out=changer.fullPivHouseholderQr().solve(model_wf);	
		if(!model_wf.isApprox(changer*out)) cout<<"bad solution!"<<endl;
		cout<<"time to solve matrix "<<(float)(clock()-t)/(CLOCKS_PER_SEC)<<endl;
		cout<<"rank: "<<lu.rank()<<" compared to: "<<head_indices.size()<<endl;
		if(!shrink &&  !post_shrink) cout<<"output symmetry: "<<(T1odd*out+(double)lil_sign(Ne)*Tphase*out).norm()<<endl;
		if(!shrink &&  !post_shrink) cout<<"output symmetry: "<<(Tinvert*out-out).norm()<<endl;
		double outnorm=sqrt(out.squaredNorm());
		double overallphase=10;
		complex<double> printout;
		for(int i=0;i<(signed)out.size();i++){
			if(norm(out(i)/outnorm)>1e-10){
				if(overallphase==10) overallphase=-arg(out(i));
				printout=out(i)/outnorm*polar(1.,overallphase);
				if(abs(printout.real())>1e-10) cout<<setprecision(14)<<printout.real()<<"\t ";
				else cout<<"0\t\t ";
				if(abs(printout.imag())>1e-10) cout<<setprecision(14)<<printout.imag()<<"\t ";
				else cout<<"0\t\t ";
				if(post_shrink) cout<<(bitset<12>)states[family_heads[i]]<<" "<<endl;
				else cout<<(bitset<12>)states[i]<<endl;
			}
		}
	}//COM for
}

//set up the change of basis matrix
void ModelTorus::makeChanger(){


}
//evaluate a single theta function for a single-particle state
complex<double> ModelTorus::landau_basis(int m){
	int type=4;
	if(NPhi%2==1) type=1;
	complex <double> tau=complex<double>(0,1./NPhi);
	complex <double> z=complex<double>(-m*M_PI/NPhi,0);
	complex<double> out;
	int sum=0;
	jacobi_theta_(&type,&z,&tau,&out,&sum);
	return out;
}

//evaluates the torus wavefunctions I am used to for comparison with duncans theta functions
double ModelTorus::my_basis(int m){
	double out=0;
	double temp;
	double eps=1e-10;
	int sign=1;
	for(int p=0; p>-1;p++){
		temp=exp(-0.5*pow(Ly*(0.5-(1.*m)/(1.*NPhi)-p),2));
		if (temp<eps) break;
		if (abs(m+p*NPhi)%2==1) sign=-1;
		else sign=1;
		out+=temp*sign;
	}
	for(int p=-1; p<0;p--){
		temp=exp(-0.5*pow(Ly*(0.5-(1.*m)/(1.*NPhi)-p),2));
		if (temp<eps) break;
		if (abs(m+p*NPhi)%2==1) sign=-1;
		else sign=1;
		out+=temp*sign;
	}
	return out*sqrt(NPhi);
}
complex<double> ModelTorus::many_body_laughlin(vector<int> sites){
	//center of mass part
	double y_total=0;
	for(int i=0;i<(signed)sites.size();i++){
		y_total+=sites[i]*Ly/(1.*NPhi);
	}
	complex<double> z,sigma;
	int reduce=1;	
	int invNu=NPhi/Ne;
//	cout<<"com part"<<endl;
	complex<double> com_part(1,0);
	for(int i=1;i<=invNu;i++){
		z=complex<double>( 0,y_total- ( (i-0.5)/(1.*invNu)-0.5)*Ly);
		weierstrass_sigma_(&z,&L1,&L2,&sigma,&reduce);
		com_part*=sigma;
//		cout<<z<<" "<<sigma<<" "<<com_part<<endl;
	}
//	com_part*=(double)(lil_sign(accumulate( sites.begin(),sites.end(), 0)/NPhi))  * exp(-invNu/(4.*NPhi)*pow(y_total,2));
	com_part*=exp(-invNu/(4.*NPhi)*pow(y_total,2));
	
	//jastrow part
//	cout<<"jastrow part ";
	complex<double> jastrow_part(1,0);
	double y_diff;
	for(int i=0;i<(signed)sites.size();i++){
		for(int j=i+1;j<(signed)sites.size();j++){
			y_diff=(sites[i]-sites[j])*Ly/(1.*NPhi);
			z=complex<double>(0,y_diff);
			weierstrass_sigma_(&z,&L1,&L2,&sigma,&reduce);
//			cout<<sigma<<endl;
//			cout<<z<<" "<<sigma<<endl;
			jastrow_part*=pow(sigma,invNu); 
			jastrow_part*=exp(-invNu*pow(y_diff,2)/(4.*NPhi));
		}
	}
//	cout<<com_part<<" "<<jastrow_part<<endl;
	return jastrow_part*com_part;
}
complex<double> ModelTorus::duncan_cfl(const vector<int> &sites){
	//center of mass part
	double y_total=0;
	for(int i=0;i<Ne;i++){
		y_total+=sites[i]*Ly/(1.*NPhi);
	}
	complex<double> z,sigma;
	int reduce=1;	
//	cout<<"com part"<<endl;
	complex<double> com_part(1,0),com_prefactor(1,0);
	complex<double> temp;
	for(int i=1;i<=invNu;i++){
		z=complex<double>( 0,y_total- ( (i-0.5)/(1.*invNu)-0.5)*Ly)-dsum/(1.*invNu);
		weierstrass_sigma_(&z,&L1,&L2,&sigma,&reduce);
//		cout<<"com_part: "<<sigma*exp(-1./(4.*NPhi)*norm(z))<<" "<<z<<" "<<Ly<<endl;
		com_prefactor*=exp(-1./(4.*NPhi)*norm(z));
		com_part*=sigma;
//		cout<<z<<" "<<sigma<<" "<<com_part<<endl;
	}
	//jastrow part
//	cout<<"jastrow part ";
	complex<double> element,jastrow_part;
	Eigen::MatrixXcd M(Ne,Ne);
	for(int i=0;i<Ne;i++){
		for(int j=0;j<Ne;j++){
			element=complex<double>(1,0);
			for(int k=0;k<Ne;k++){
				if(k==i) continue;
				z=complex<double>(0,(sites[i]-sites[k])*Ly/(1.*NPhi)) - ds[j]+dsum/(1.*Ne);
				weierstrass_sigma_(&z,&L1,&L2,&sigma,&reduce);
				element*=sigma;
			}
			//M(i,j)=element*exp(-sites[i]*Ly*ds[i].imag()/(2.*invNu*NPhi));
			M(i,j)=element*exp( ((double)(sites[i])*L2/(1.*NPhi) ) *conj(ds[j])/(2.*invNu));
			//M(i,j)=element*exp( ((double)(sites[i])*conj(L2)/(1.*NPhi) ) *ds[j]/(2.*invNu));
		}
	}
//	cout<<M<<endl;
	Eigen::FullPivLU<Eigen::MatrixXcd> Msolver(M);
	//compute the (complicated) exponential prefactor to the determinant
	complex<double> sum(0,0);
	double sum2=0;
	complex<double> dvar(0,0);
	for(int i=0;i<Ne;i++){
		sum+=0.25*pow( sites[i]*Ly/(1.*NPhi) , 2);
		dvar+=norm(dsum/(1.*Ne)-ds[i]);
		for(int k=0;k<Ne;k++){
			if(k==i) continue;
			sum2+=pow( (sites[i]-sites[k])*Ly/(1.*NPhi),2);
		}
	}
	jastrow_part=Msolver.determinant();

	//cout<<com_part*com_prefactor<<" "<<jastrow_part*exp(-sum2/(4.*NPhi))*exp(-(double)(Ne-1)*dvar/(4.*NPhi))*exp(-1./(2.*NPhi)*complex<double>(0,y_total)*dsum)<<endl;
	return jastrow_part*com_part*exp(-sum);
}
complex<double> ModelTorus::single_determinant(int k, const vector<int> &eval_sites_element, const vector< complex<double> > &single_particle_orbitals){
	Eigen::FullPivLU<Eigen::MatrixXcd> LUsolver;	
	vector<int> eval_states_element=bitset_to_pos(states[k],NPhi);	
	Eigen::MatrixXcd to_determinant=Eigen::MatrixXcd::Zero(Ne,Ne);

	for(int n=0;n<Ne;n++){
		for(int m=0;m<Ne;m++){
			if(eval_sites_element[n]>eval_states_element[m]) to_determinant(n,m)=(double)(lil_sign(eval_states_element[m]))*single_particle_orbitals[eval_sites_element[n]-eval_states_element[m]];
			else to_determinant(n,m)=(double)(lil_sign(eval_states_element[m])*lil_sign(NPhi))*single_particle_orbitals[-eval_sites_element[n]+eval_states_element[m]];
		}
	}
	//compute determinant
//	cout<<to_determinant<<endl<<endl;
	LUsolver.compute(to_determinant);
	return LUsolver.determinant();
}
complex<double> ModelTorus::duncan_test(int m){
	complex<double> z(0,m*Ly/NPhi),sigma;
	int reduce=1;
	weierstrass_sigma_(&z,&L1,&L2,&sigma,&reduce);
	return sigma*exp(-0.25*pow(m*Ly/NPhi,2)/(1.*NPhi));
}
void ModelTorus::make_translate_one(){
	T1odd=Eigen::Matrix<double,-1,-1>::Zero(nStates,nStates);
	T1even=Eigen::Matrix<double,-1,-1>::Zero(nStates,nStates);
	Tinvert=Eigen::Matrix<double,-1,-1>::Zero(nStates,nStates);
	vector<int> sites;
	int newstate=0,newindex;
	vector<int>::const_iterator it;
	for(int i=0;i<nStates;i++){

		//translation matrix
		newstate=0;
		sites=bitset_to_pos(states[i],NPhi);
		for(int j=0;j<(signed)sites.size(); j++) newstate=newstate | 1<<((sites[j]+1)%NPhi);
		it=find(states.begin(),states.end(),newstate);
		if(it!=states.end()){
			newindex=it-states.begin();
			T1even(i,newindex)=1;
			T1odd(i,newindex)=1;
			if( (states[newindex] <= states[i]) && Ne%2==0){
//					cout<<"wrapped states "<<(bitset<9>)states[newindex]<<" "<<(bitset<9>)states[i]<<endl;
				T1odd(i,newindex)*=-1;
			}
		}else{
			cout<<"state: "<<(bitset<9>)newstate<<" not found"<<endl;
			exit(0);
		}
		
		//inversion matrix
		newstate=0;
		sites=bitset_to_pos(states[i],NPhi);
		for(int j=0;j<(signed)sites.size();j++){
			newstate=newstate | 1<<( (NPhi-1-sites[j])%NPhi);
		}
		it=find(states.begin(),states.end(),newstate);
		if(it!=states.end()){
			newindex=it-states.begin();
			Tinvert(i,newindex)=1;
		}else{
			cout<<"state: "<<(bitset<9>)newstate<<" not found"<<endl;
			exit(0);
		}
		
	}
}
void ModelTorus::make_states(){
	int compare_bit;
	bool found;
	vector<int>::iterator it;
	for(int i=0;i<1<<NPhi;i++){
		if(count_bits(i)==Ne){
			if(!shrink) states.push_back(i);
			else{
				compare_bit=i;
				found=false;
//				cout<<"searching for"<<(bitset<6>)i<<endl;
				for(int j=0;j<NPhi-1;j++){
//					cout<<(bitset<6>)compare_bit<<endl;
					compare_bit=cycle_bits(compare_bit,NPhi);
					it=find(states.begin(),states.end(),compare_bit);
					if(it!=states.end()){
//						cout<<"found translated!"<<endl;
						family_size[it-states.begin()]++;
						found=true;
						break;
					}
				}
				if(!found){
					compare_bit=invert_bits(i,NPhi);
					for(int j=0;j<NPhi;j++){
						it=find(states.begin(),states.end(),compare_bit);
						if(it!=states.end()){
	//						cout<<"found translated!"<<endl;
							family_size[it-states.begin()]++;
							invert_family[it-states.begin()]=true;
							break;
						}
						compare_bit=cycle_bits(compare_bit,NPhi);
						if(j==NPhi-1){
							states.push_back(i);
							family_size.push_back(1);
							invert_family.push_back(false);
						}					
					}//for through shifts
				}				
			}//else shrink
		}
	}
	nStates=states.size();
//	cout<<"nStates: "<<nStates<<endl;
//	for(int i=0;i<nStates;i++){
//		cout<<(bitset<12>)states[i];
//		if(shrink) cout<<" "<<family_size[i]<<" "<<invert_family[i];
//		cout<<endl;
//	}
	
}
void ModelTorus::translation_eigenvectors(){
	vector <int> heads;
	vector < vector<int> > families;
	vector<int> tempvec;
	int tempstate;
	bool found;
	for(int i=0;i<nStates;i++){
		found=false;
		cout<<"searching for"<<(bitset<8>)states[i]<<endl;
		//if this state is already in a family, add it to that family
		for(int j=0;j<(signed)families.size();j++){
			tempstate=heads[j];
			for(int k=0;k<NPhi-1;k++){
				tempstate=cycle_bits(tempstate,NPhi);
				//cout<<"comparing to"<<(bitset<8>)tempstate<<endl;
				if(tempstate==states[i]){
					families[j].push_back(states[i]);
					cout<<"match found"<<endl;
					found=true;
					break;
				}
			}
			if(found) break;
		}
		//it not, make a new family
		if(!found){
			cout<<"its a new family"<<endl;
			tempvec.clear();
			tempvec.push_back(states[i]);
			families.push_back(tempvec);
			heads.push_back(states[i]);
		}		
	}
	for(int i=0;i<(signed)families.size();i++){
		for(int j=0;j<(signed)families[i].size();j++){
			cout<<(bitset<8>)families[i][j]<<endl;
		}
		cout<<endl;
	} 
}
//takes a bit and a set of translation family heads, and finds the translation family the bit belongs to

int ModelTorus::find_translated(int bit, const vector<int> &heads){
	int sign=1;
	int teststate,newteststate;
	for (int i=0; i<(signed)heads.size(); i++){
		sign=1;
		teststate=states[heads[i]];
		for(int j=0; j<NPhi; j++){
			if(bit==teststate) return i*sign;
			newteststate=cycle_bits(teststate,NPhi);
			sign*=-1;
			if(newteststate<teststate) sign*=-1;
			teststate=newteststate;
		}
	}
	return -1000;
}	
int main(){

	ModelTorus MT;
	MT.basis_change();
}	

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "MOint.h"
#include "permute.h"
#include "slaterd.h"

/*Matrix Decomposition Full Configuration Interaction
 Based on Koch, Henrik, and Esper Dalgaard. "A variational matrix decomposition applied to full configuration-interaction calculations." Chemical physics letters 198.1 (1992): 51-58.

Most of the functionality is an implementation of equation 24 and its counterpart in beta
*/



INIT_PLUGIN

using namespace boost;

namespace psi{ namespace variational_decomp_fci {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "VARIATIONAL_DECOMP_FCI"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

void normalize(double *vec, int length);
double randouble();
void test_factored_sigma(int n_terms,int alphae, int betae, int nmo, double eSCF,double **mo_OEIprime, double *mo_TEI);
void variational_matrix_decomposition(int state, int nmo, int alphae, int betae, double eSCF, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,int print);
int davidP(int state, int M, double *total_energy,double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,int use_guess,int print);
int davidQ(int state, int M, double *total_energy,double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,int use_guess,int print);
double get_energyP(int state, double **P, double **sigmaP, int n_Aterms, int alphae, int nmo);
double get_energyQ(int state, double **Q, double **sigmaQ, int n_Bterms, int betae, int nmo);
void get_factored_sigmaP(double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print);
void get_factored_sigmaQ(double *P,double *Q,double *Pprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print);
void get_factored_sigmaPAA(double *P,double *Q,double *Qprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaPBB(double *P,double *Q,double *Qprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaPAB(double *P,double *Q,double *Qprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQAA(double *P,double *Q,double *Pprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQBB(double *P,double *Q,double *Pprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_factored_sigmaQAB(double *P,double *Q,double *Pprime,double *sigmaP,  int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI);
void get_augmented_sigmaP(double *Pnew, double *Qnew, double **P, double **Q, double *sigmaP, double E0, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI, int print);
void get_augmented_sigmaQ(double *Pnew, double *Qnew, double **P, double **Q, double *sigmaP, double E0, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI, int print);
void test_augmented_sigma(int n_terms, int alphae, int betae, int nmo, double eSCF, double **mo_OEIprime, double *mo_TEI);
extern "C" 
PsiReturnType variational_decomp_fci(Options& options)
{
    int print = options.get_int("PRINT");
/* INITIALIZATIONS */
	//pointers to important constants
	int *nmo_ = new int; //# molecular orbitals
	int *alphae_ = new int; //# alpha electrons
	int *betae_ = new int; //# beta electrons
	double *nuc_rep_energy_= new double;
	double *eSCF_ = new double;
	const char *molname = new char[80], *basisname = new char[80]; //molecule name and basis name

	//find important constants
	//function in MOint.h
	initialize_MO_constants(nmo_,alphae_,betae_,nuc_rep_energy_,eSCF_,molname,basisname);
	//recast pointers as ints
	int nmo = *(nmo_);
	int alphae = *alphae_;
	int betae = *betae_;
	double nuc_rep_energy= *(nuc_rep_energy_);
	double eSCF = *(eSCF_);
	eSCF-=nuc_rep_energy;
	delete [] eSCF_;
	delete [] nmo_;
	delete [] alphae_;
	delete [] betae_;

	//initialize one-electron and two-electron integrals
	double **mo_OEI=block_matrix(nmo,nmo); //Pointer to OEI in MO basis
	double *mo_TEI=init_array(nmo*nmo*nmo*nmo);

	//transform MOs, storing in **OEI and *TEI arrays
	//function in MOint.h
	MO_transform(mo_OEI,mo_TEI,nmo);

	//number of determinants
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);                                                          
	//full size of one dimension of NxN Hamiltonian  
	int N = astringcount*bstringcount;     

	//option to build and diagonalize the full Hamiltonian to test convergence code
	double etot; //store ground state energy
	int buildH = 1;
        if(buildH){                       
		double **H = block_matrix(astringcount*bstringcount,astringcount*bstringcount);
		build_full_Hamiltonian(H,mo_OEI,mo_TEI,alphae,betae,nmo);
		double *lambda = init_array(astringcount*bstringcount);
		double *work= init_array(10*astringcount*bstringcount);
		C_DSYEV('V','U',astringcount*bstringcount,&(H[0][0]),astringcount*bstringcount,lambda,work,10*astringcount*bstringcount);
		for(int i=0;i<astringcount*bstringcount;i++){
			fprintf(outfile,"%lf\n",lambda[i]);
		}
		etot = lambda[0];
		free(lambda);
		free(work);
		free_block(H);
	}
	
	//copy OEI to new array to be updated with TEI terms
	double **mo_OEIprime = block_matrix(nmo,nmo);
	C_DCOPY(nmo*nmo,&(mo_OEI[0][0]),1,&(mo_OEIprime[0][0]),1);	

	//add TEI to OEI to build matrix of hkl'
	for(int k =0;k<nmo;k++){
		for(int l=0;l<nmo;l++){
			for(int j =0;j<nmo;j++){
				mo_OEIprime[k][l]-=0.5*mo_TEI[k*nmo*nmo*nmo+j*nmo*nmo+j*nmo+l];
			}
		}
	}
	print_mat(mo_OEIprime,nmo,nmo,outfile);

	//functions to test sigma builds and augmented sigma builds with the c0 normalization parameter	
	int test_sig = 1;
	if(test_sig){
		//test_factored_sigma(1,alphae, betae, nmo, eSCF,mo_OEIprime, mo_TEI);
		test_augmented_sigma(2,  alphae, betae,  nmo, eSCF, mo_OEIprime, mo_TEI);
	}

	//back-and-forth Davidson
	//variational_matrix_decomposition(0,nmo,alphae,betae,eSCF,mo_OEI,mo_OEIprime,mo_TEI,2);

        if(buildH){       
		printf("FCI Energy = %lf\n",etot);
	}                
	free(mo_TEI);
	free_block(mo_OEIprime);
	
/* END INITIALIZATIONS */
return Success;
}

//return random double
double randouble(){
	double F = (double)rand() / RAND_MAX;
	return F;
}

//tests sigma builds augmented with c0 normalization parameter (eq. 18)
void test_augmented_sigma(int n_terms, int alphae, int betae, int nmo, double eSCF, double **mo_OEIprime, double *mo_TEI){
	printf("\n--- test_augmented_sigma ---\n");
	fprintf(outfile,"+++++++++++++++++++++++ TEST AUGMENTED SIGMA ++++++++++++++++++++++++++++++++\n");
	int i,j,k,l;
	//number of P and Q tables
	n_terms = 2;

	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int NA = astringcount+1;
	int NB = bstringcount+1;
	//number of guesses for current, un-converged iteration
	int n_ortho = 4;

	//P tables are indexed so that P[0] = c0, the normalization parameter
	//Q are the same
	double **P = block_matrix(n_ortho,NA);
	double **sigmaP = block_matrix(n_ortho,NA);
	double **Q = block_matrix(n_ortho,NB);
	double **sigmaQ = block_matrix(n_ortho,NB);

	double **ortho_quadP = block_matrix(n_ortho,n_ortho);
	double **ortho_quadQ = block_matrix(n_ortho,n_ortho);
	double **overlapP = block_matrix(n_ortho,n_ortho);
	double **overlapQ = block_matrix(n_ortho,n_ortho);
	
	srand(time(NULL));
	
	//set "converged" vector to Hartree-Fock guess
	P[0][1] = 1.0;
	//populate guess P tables with random values schmidt orthogonalized against converged vector and each other
	for(int i=1;i<n_ortho;i++){
		for(int j=0;j<NA;j++){
			P[i][j]=randouble();
		}
		normalize(&(P[i][0]),NA);
		for(int j=0;j<i;j++){
			double proj = C_DDOT(NA,&(P[i][0]),1,&(P[j][0]),1);
			C_DAXPY(NA,-proj,&(P[j][0]),1,&(P[i][0]),1);
		}
		normalize(&(P[i][0]),NA);
	}
	//if the system is a singlet, copy P over to Q tables
	//otherwise, choose new random vectors orthogonalized against each other
	if(astringcount == bstringcount){
		for(i=0;i<n_ortho;i++){
			for(j=0;j<NB;j++){
				Q[i][j] = P[i][j];
			}
		}
	}
	else{
		Q[0][1] = 1.0;
		for(int i=1;i<n_ortho;i++){
			for(int j=0;j<NB;j++){
				Q[i][j]=randouble();
			}
			normalize(&(Q[i][0]),NB);
			for(int j=0;j<i;j++){
				double proj = C_DDOT(NB,&(Q[i][0]),1,&(Q[j][0]),1);
				C_DAXPY(NB,-proj,&(Q[j][0]),1,&(Q[i][0]),1);
			}
			normalize(&(Q[i][0]),NB);
		}
	}
	//compute augmented sigma for all guess vectors wrt to converged vector
	for(i=1;i<n_ortho;i++){
		get_augmented_sigmaP(&(P[i][0]), &(Q[i][0]), P, Q, &(sigmaP[i][0]), eSCF, n_terms, n_terms, alphae, betae,  nmo, mo_OEIprime, mo_TEI, 1);
	}
	//compute overlap of guess vectors	
	C_DGEMM('n','t', n_ortho, n_ortho, NA, 1.0, &(P[0][0]), NA,&(P[0][0]), NA, 0.0, &(overlapP[0][0]), n_ortho);
	//compute quadratic form P*sigma^P
	C_DGEMM('n','t', n_ortho, n_ortho, NA, 1.0, &(P[0][0]), NA,&(sigmaP[0][0]), NA, 0.0, &(ortho_quadP[0][0]), n_ortho);
	
	//naive checks:
	//check overlap of vectors
	double ovP_orth = 0.0;
	int ovP_norm = 0;
	for(i=0;i<n_ortho;i++){
		if(fabs(overlapP[i][i] - 1.0) > 10E-6){
			ovP_norm = 1;
		}
		for(j=0;j<i;j++){
			ovP_orth += overlapP[i][j] - overlapP[j][i];
		}
	}
	if(!ovP_norm){
		printf("P tables ARE normalized\n");
	}
	else{
		printf("P tables ARE NOT normalized\n");
	}
	if(fabs(ovP_orth) < 10E-6){
		printf("P tables ARE orthogonal\n");
	}
	else{
		printf("P table ARE NOT orthogonal\n");
	}
	
	double sigP_symm = 0.0;
	for(i=1;i<n_ortho;i++){
		for(j=1;j<i;j++){
			sigP_symm += ortho_quadP[i][j] - ortho_quadP[j][i];
		}
	}
	if(fabs(sigP_symm) < 10E-6){
		printf("P quadratic form IS symmetric \n");
	}
	else{
		printf("P quadratic form IS NOT symmetric: TEST FAILED \n");
	}
	fprintf(outfile,"		/// P ///\n");
	fprintf(outfile,"test P %d\n",n_ortho);
	print_mat(P,n_ortho,NA,outfile);
	fprintf(outfile,"overlap P %d\n",n_ortho);
	print_mat(overlapP,n_ortho,n_ortho,outfile);
	fprintf(outfile,"test sigmaP %d\n",n_ortho);
	print_mat(sigmaP,n_ortho,NA,outfile);
	fprintf(outfile,"P quad %d\n",n_ortho);
	print_mat(ortho_quadP,n_ortho-1,n_ortho-1,outfile);

	
	//compute augmented sigma^Q
	for(i=1;i<n_ortho;i++){
		get_augmented_sigmaQ(&(P[i][0]), &(Q[i][0]), P, Q, &(sigmaQ[i][0]), eSCF, n_terms, n_terms, alphae, betae,  nmo, mo_OEIprime, mo_TEI, 1);
	}
	
	C_DGEMM('n','t', n_ortho, n_ortho, NB, 1.0, &(Q[0][0]), NB,&(Q[0][0]), NB, 0.0, &(overlapQ[0][0]), n_ortho);
	C_DGEMM('n','t', n_ortho, n_ortho, NB, 1.0, &(Q[0][0]), NB,&(sigmaQ[0][0]), NB, 0.0, &(ortho_quadQ[0][0]), n_ortho);
	double ovQ_orth = 0.0;
	int ovQ_norm = 0;
	for(i=0;i<n_ortho;i++){
		if(fabs(overlapQ[i][i] - 1.0) > 10E-6){
			ovQ_norm = 1;
		}
		for(j=0;j<i;j++){
			ovQ_orth += overlapQ[i][j] - overlapQ[j][i];
		}
	}
	if(!ovQ_norm){
		printf("Q tables ARE normalized\n");
	}
	else{
		printf("Q tables ARE NOT normalized\n");
	}
	if(fabs(ovQ_orth) < 10E-6){
		printf("Q tables ARE orthogonal\n");
	}
	else{
		printf("Q table ARE NOT orthogonal\n");
	}
	
	double sigQ_symm = 0.0;
	for(i=1;i<n_ortho;i++){
		for(j=1;j<i;j++){
			sigQ_symm += ortho_quadQ[i][j] - ortho_quadQ[j][i];
		}
	}
	if(fabs(sigQ_symm) < 10E-6){
		printf("Q quadratic form IS symmetric \n");
	}
	else{
		printf("Q quadratic form IS NOT symmetric: TEST FAILED \n");
	}
	fprintf(outfile,"		/// Q ///\n");
	fprintf(outfile,"test Q %d\n",n_ortho);
	print_mat(Q,n_ortho,NB,outfile);
	fprintf(outfile,"overlap Q %d\n",n_ortho);
	print_mat(overlapQ,n_ortho,n_ortho,outfile);
	fprintf(outfile,"test sigmaQ %d\n",n_ortho);
	print_mat(sigmaQ,n_ortho,NB,outfile);
	fprintf(outfile,"Q quad %d\n",n_ortho);
	print_mat(ortho_quadQ,n_ortho-1,n_ortho-1,outfile);
	
	free_block(overlapP);
	free_block(overlapQ);
	free_block(P);
	free_block(Q);
	free_block(sigmaP);
	free_block(sigmaQ);
	free_block(ortho_quadP);
	free_block(ortho_quadQ);
	
	printf("------\n");
}
void test_factored_sigma(int n_terms,int alphae, int betae, int nmo, double eSCF,double **mo_OEIprime, double *mo_TEI){
	int i,j,k,l;
	printf("--- test_factored_sigma ---\n");
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int NA = astringcount+1;
	int NB = bstringcount+1;
	//number of guess vectors
	int n_ortho = 2;

	double **P = block_matrix(n_ortho,NA);
	double **sigmaP = block_matrix(n_ortho,NA);
	double **Q = block_matrix(n_ortho,NB);
	double **sigmaQ = block_matrix(n_ortho,NB);

	double **ortho_quadP = block_matrix(n_ortho,n_ortho);
	double **ortho_quadQ = block_matrix(n_ortho,n_ortho);
	double **overlapP = block_matrix(n_ortho,n_ortho);
	double **overlapQ = block_matrix(n_ortho,n_ortho);
	//populate P randomly, orthonormalizing subsequent vectors
		srand(time(NULL));
		for(int i=1;i<NA;i++){
			P[0][i] = randouble();
		}
		//P[0][1]=1.0;
		normalize(&(P[0][0]),NA);
		for(int i=1;i<n_ortho;i++){
			for(int j=1;j<NA;j++){
				P[i][j]=randouble();
			}
			normalize(&(P[i][0]),NA);
			for(int j=0;j<i;j++){
				double proj = C_DDOT(NA-1,&(P[i][1]),1,&(P[j][1]),1);
				C_DAXPY(NA-1,-proj,&(P[j][1]),1,&(P[i][1]),1);
			}
			normalize(&(P[i][0]),NA);
		}
	//if a singlet, copy P to Q
	//otherwise, repeat ON process
		if(astringcount == bstringcount){
			for(int i=0;i<n_ortho;i++){
				for(int j=0;j<NA;j++){
					Q[i][j]=P[i][j];
				}
			}
		}
		else{
			Q[0][1]=1.0;
			normalize(&(Q[0][0]),NB);
			for(int i=1;i<n_ortho;i++){
				for(int j=1;j<NB;j++){
					Q[i][j]=randouble();
				}
				normalize(&(Q[i][0]),NB);
				for(int j=0;j<i;j++){
					double proj = C_DDOT(NB-1,&(Q[i][1]),1,&(Q[j][1]),1);
					C_DAXPY(NB-1,-proj,&(Q[j][1]),1,&(Q[i][1]),1);
				}
				normalize(&(Q[i][0]),NB);
			}
		}

		fprintf(outfile,"+++++++++++++++++++++++ TEST FACTORED SIGMA ++++++++++++++++++++++++++++++++\n");
		for(int i=0;i<n_ortho;i++){
			get_factored_sigmaP(&(P[i][1]),&(Q[1][1]),&(Q[1][1]),&(sigmaP[i][1]), alphae, betae, nmo, mo_OEIprime, mo_TEI,3);
		}
		C_DGEMM('n','t', n_ortho, n_ortho, NA, 1.0, &(P[0][0]), NA,&(P[0][0]), NA, 0.0, &(overlapP[0][0]), n_ortho);
		C_DGEMM('n','t', n_ortho, n_ortho, NA, 1.0, &(P[0][0]), NA,&(sigmaP[0][0]), NA, 0.0, &(ortho_quadP[0][0]), n_ortho);
		
		//naive checks:
		//check overlap of vectors
		double ovP_orth = 0.0;
		int ovP_norm = 0;
		printf("-- check P:\n");
		for(i=0;i<n_ortho;i++){
			if(fabs(overlapP[i][i] -  1.0) > 10E-6){
				ovP_norm = 1;
			}
			for(j=0;j<i;j++){
				ovP_orth += overlapP[i][j] - overlapP[j][i];
			}
		}
		if(!ovP_norm){
			printf("P tables ARE normalized\n");
		}
		else{
			printf("P tables ARE NOT normalized\n");
		}
		if(fabs(ovP_orth) < 10E-6){
			printf("P tables ARE orthogonal\n");
		}
		else{
			printf("P table ARE NOT orthogonal\n");
		}
		
		double sigP_symm = 0.0;
		for(i=1;i<n_ortho;i++){
			for(j=1;j<i;j++){
				sigP_symm += ortho_quadP[i][j] - ortho_quadP[j][i];
			}
		}
		if(fabs(sigP_symm) < 10E-6){
			printf("P quadratic form IS symmetric \n");
		}
		else{
			printf("P quadratic form IS NOT symmetric: TEST FAILED \n");
		}
		printf("--\n");
		fprintf(outfile,"		/// P ///\n");
		fprintf(outfile,"test P %d\n",n_ortho);
		print_mat(P,n_ortho,NA,outfile);
		fprintf(outfile,"overlap P %d\n",n_ortho);
		print_mat(overlapP,n_ortho,n_ortho,outfile);
		fprintf(outfile,"test sigmaP %d\n",n_ortho);
		print_mat(sigmaP,n_ortho,NA,outfile);
		fprintf(outfile,"P quad %d\n",n_ortho);
		print_mat(ortho_quadP,n_ortho,n_ortho,outfile);
		double *lambda = init_array(n_ortho);
		double *work= init_array(10*n_ortho);
		C_DSYEV('V','U',n_ortho,&(ortho_quadP[0][0]),n_ortho,lambda,work,10*n_ortho);
		//print_mat(ortho_quadP,n_ortho,n_ortho,outfile);
		for(int i=0;i<n_ortho;i++){
			fprintf(outfile,"%lf\n",lambda[i]);
		}
		for(int i=0;i<n_ortho;i++){
			get_factored_sigmaQ(&(P[1][1]),&(Q[i][1]),&(P[1][1]),&(sigmaQ[i][1]), alphae, betae, nmo, mo_OEIprime, mo_TEI,3);
		}
		C_DGEMM('n','t', n_ortho, n_ortho, NB, 1.0, &(Q[0][0]), NB,&(Q[0][0]), NB, 0.0, &(overlapQ[0][0]), n_ortho);
		C_DGEMM('n','t', n_ortho, n_ortho, NB, 1.0, &(Q[0][0]), NB,&(sigmaQ[0][0]), NB, 0.0, &(ortho_quadQ[0][0]), n_ortho);
	double ovQ_orth = 0.0;
	int ovQ_norm = 0.0;
	printf("-- check Q:\n");
	for(i=0;i<n_ortho;i++){
		if(fabs(overlapQ[i][i] - 1.0) > 10E-6){
			ovQ_norm = 1;
		}
		for(j=0;j<i;j++){
			ovQ_orth += overlapQ[i][j] - overlapQ[j][i];
		}
	}
	if(!ovQ_norm){
		printf("Q tables ARE normalized\n");
	}
	else{
		printf("Q tables ARE NOT normalized\n");
	}
	if(fabs(ovQ_orth) < 10E-6){
		printf("Q tables ARE orthogonal\n");
	}
	else{
		printf("Q table ARE NOT orthogonal\n");
	}
	
	double sigQ_symm = 0.0;
	for(i=1;i<n_ortho;i++){
		for(j=1;j<i;j++){
			sigQ_symm += ortho_quadQ[i][j] - ortho_quadQ[j][i];
		}
	}
	if(fabs(sigQ_symm) < 10E-6){
		printf("Q quadratic form IS symmetric \n");
	}
	else{
		printf("Q quadratic form IS NOT symmetric: TEST FAILED \n");
	}
	printf("--\n");


		fprintf(outfile,"		/// Q ///\n");
		fprintf(outfile,"test Q %d\n",n_ortho);
		print_mat(Q,n_ortho,NB,outfile);
		fprintf(outfile,"overlap Q %d\n",n_ortho);
		print_mat(overlapQ,n_ortho,n_ortho,outfile);
		fprintf(outfile,"test sigmaQ %d\n",n_ortho);
		print_mat(sigmaQ,n_ortho,NB,outfile);
		fprintf(outfile,"Q quad %d\n",n_ortho);
		print_mat(ortho_quadQ,n_ortho,n_ortho,outfile);
		free_block(overlapP);
		free_block(overlapQ);
		free_block(P);
		free_block(Q);
		free_block(sigmaP);
		free_block(sigmaQ);
		free_block(ortho_quadP);
		free_block(ortho_quadQ);
}


//runs variational matrix decomposition algorithm by adding one alpha and one beta determinant at a time
void variational_matrix_decomposition(int state, int nmo, int alphae, int betae, double eSCF, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,int print){

	int i,j,k,l;

	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	
	int wfn_terms, dav_iters;
	int max_wfn_terms = 1;
	int max_dav_iters = 100;

	//number of Davidson roots
	int roots = 1;

	//converged tables
	double **P = block_matrix(roots,max_wfn_terms*astringcount);
	double **Q = block_matrix(roots,max_wfn_terms*bstringcount);

	//energy tolerance for convergence
	double Etol = 10E-7;

	//best guesses for energy from P and Q iterations
	double *Ep = init_array(roots);
	double *Eq = init_array(roots);
	//initialize total energy to HF energy
	Ep[0]=eSCF;
	Eq[0]=eSCF;

	//current iteration energies
	double energyP = 1.0;	
	double energyQ = 0.0;

	//converged energies from all iterations	
	double *conv_energies = init_array(max_wfn_terms);

	for(wfn_terms=0;wfn_terms<max_wfn_terms;wfn_terms++){
		for(i=0;i<astringcount;i++){
			P[state][wfn_terms*astringcount+i]=1/sqrt(astringcount);
		}
		for(i=0;i<wfn_terms;i++){
			double proj = C_DDOT(astringcount,&(P[state][wfn_terms*astringcount]),1,&(P[state][i*astringcount]),1);
			C_DAXPY(astringcount,-proj,&(P[state][i*astringcount]),1,&(P[state][wfn_terms*astringcount]),1);
		}
		normalize(&(P[state][wfn_terms*astringcount]),astringcount);
		//need to guess a table for Q, so choose homogeneous guess schmidt orthogonalized against the previous vectors
		for(i=0;i<bstringcount;i++){
			Q[state][wfn_terms*bstringcount+i]=1/sqrt(bstringcount);
		}
		for(i=0;i<wfn_terms;i++){
			double proj = C_DDOT(bstringcount,&(Q[state][wfn_terms*bstringcount]),1,&(Q[state][i*bstringcount]),1);
			C_DAXPY(bstringcount,-proj,&(Q[state][i*bstringcount]),1,&(Q[state][wfn_terms*bstringcount]),1);
		}
		normalize(&(Q[state][wfn_terms*bstringcount]),bstringcount);
		dav_iters = 0;	
		fprintf(outfile,"---------------------------WFN TERMS = %d---------------------------\n",wfn_terms);
		//this is a switch for using a guess internally in the Davidson function (use_guess = 0) or using the guess from the previous Davidson iteration
		//typically, use the internal guess for the first iteration and reuse the previous iteration's result for subsequent iterations
		int use_guess = 0;
		//iterate over P and Q, holding each fixed in turn until their energies are self consistent
		while(fabs(energyP-energyQ)>Etol && dav_iters < max_dav_iters){
			//print current tables
			fprintf(outfile,"P tables\n");
			print_mat(P,1,(wfn_terms+1)*astringcount,outfile);
			fprintf(outfile,"Q tables\n");
			print_mat(Q,1,(wfn_terms+1)*bstringcount,outfile);
			printf("----------------------------------------------------------\n");
			printf("P DAVIDSON ITERATION = %d\n",dav_iters);
			//Q is fixed
			davidP(0,roots, Ep,P, Q, wfn_terms+1, wfn_terms+1, alphae, betae, nmo, mo_OEI,mo_OEIprime, mo_TEI,use_guess,1);
			//store ground state energy
			energyP = Ep[0];
			//print P table
			if(print>1){
				printf("P wavefunction\n");
				for(j=0;j<astringcount;j++){
					printf("%lf ",P[state][wfn_terms*astringcount+j]);
				}
				printf("\n");
			}
			printf("Q DAVIDSON ITERATION = %d\n",dav_iters);
			//print current tables
			fprintf(outfile,"P tables\n");
			print_mat(P,1,(i+1)*astringcount,outfile);
			fprintf(outfile,"Q tables\n");
			print_mat(Q,1,(i+1)*bstringcount,outfile);
			//P is fixed
			davidQ(0,roots, Eq,P, Q, wfn_terms+1, wfn_terms+1, alphae, betae, nmo,mo_OEI, mo_OEIprime, mo_TEI,use_guess,1);
			//store ground state energy
			energyQ = Eq[0];
			if(print>1){
				printf("Q wavefunction\n");
				for(j=0;j<bstringcount;j++){
					printf("%lf ",Q[state][wfn_terms*bstringcount+j]);
				}
				printf("\n");
			}
			dav_iters++;
			use_guess = 1;
			printf("deltaE = %lf\n",fabs(energyP-energyQ));
		}
		conv_energies[wfn_terms] = Ep[0];
		if(dav_iters<max_dav_iters){
			printf("|||| MATRIX DECOMP ITERATION %d CONVERGED IN %d DAVIDSON ITERATIONS ||||\n",wfn_terms,dav_iters);
		}
		else{
			printf("++++ MATRIX DECOMP ITERATION %d DID NOT CONVERGE ++++\n",wfn_terms);
		}
		printf("Ep = %lf\nEq = %lf\n",Ep[0],Eq[0]);
		fprintf(outfile,"P %d\n",wfn_terms+1);
		print_mat(P,1,(wfn_terms+1)*astringcount,outfile);
		fprintf(outfile,"Q %d\n",wfn_terms+1);
		print_mat(Q,1,(wfn_terms+1)*bstringcount,outfile);
		fprintf(outfile,"P norm = %20.14lf\n",C_DDOT((wfn_terms+1)*astringcount,&(P[state][0]),1,&(P[state][0]),1));
		fprintf(outfile,"Q norm = %20.14lf\n",C_DDOT((wfn_terms+1)*bstringcount,&(Q[state][0]),1,&(Q[state][0]),1));
	}
	printf("\nVARIATIONAL MATRIX DECOMPOSITION RESULTS\n");
	printf("Terms | Energy\n");
	printf("______________\n");
	for(i=0;i<max_wfn_terms;i++){
		printf("%5d | %lf\n",i,conv_energies[i]);
	}
	printf("\n");
	free(conv_energies);

}


//normalize a vector 
//if the norm is too small, the vector is set to zero
void normalize(double *vec, int length){
	double norm = C_DDOT(length,vec,1,vec,1);
	norm = sqrt(norm);
	if(norm>10E-6){
		norm = 1/norm;
		C_DSCAL(length,norm,vec,1);
	}
	else{
		C_DSCAL(length,0,vec,1);
	}
}

//compute augmented sigma^P (including first term for previous iteration's energy)
void get_augmented_sigmaP(double *Pnew, double *Qnew, double **P, double **Q, double *sigmaP, double E0, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI, int print){
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	sigmaP[0] = Pnew[0]*E0;
	for(i=0;i<n_Aterms-1;i++){
		get_factored_sigmaP(&(P[0][i*astringcount]),&(Q[0][i*bstringcount]),&Qnew[1],&sigmaP[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	}
	sigmaP[0] += C_DDOT(astringcount,&sigmaP[1],1,&Pnew[1],1);
	C_DSCAL(astringcount,Pnew[0],&sigmaP[1],1);
	get_factored_sigmaPBB(&Pnew[1],&Qnew[1],&Qnew[1],&sigmaP[1],alphae,betae,nmo,mo_OEIprime,mo_TEI);
}

//compute augmented sigma^Q (including first term for previous iteration's energy)
void get_augmented_sigmaQ(double *Pnew, double *Qnew, double **P, double **Q, double *sigmaQ, double E0, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI, int print){
	int i,j,k,l;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	sigmaQ[0] = Qnew[0]*E0;
	for(i=0;i<n_Bterms-1;i++){
		get_factored_sigmaQ(&(P[0][i*astringcount]),&(Q[0][i*bstringcount]),&Pnew[1],&sigmaQ[1],alphae,betae,nmo,mo_OEIprime,mo_TEI,0);
	}
	sigmaQ[0] += C_DDOT(bstringcount,&sigmaQ[1],1,&Qnew[1],1);
	C_DSCAL(bstringcount,Qnew[0],&sigmaQ[1],1);
	get_factored_sigmaQAA(&Pnew[1],&Qnew[1],&Pnew[1],&sigmaQ[1],alphae,betae,nmo,mo_OEIprime,mo_TEI);
}

//compute sigma^P for a given state up to a given number of alpha and beta terms in the trial wavefunction
void get_factored_sigmaP(double *P,double *Q,double *Qprime,double *sigmaP, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print){
	int astringcount = nchoosek(nmo,alphae);
	get_factored_sigmaPAA(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>0){
		printf("sigmaP 1:     ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}
	get_factored_sigmaPBB(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>1){
		printf("sigmaP 1+2:   ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}
	get_factored_sigmaPAB(P,Q,Qprime,sigmaP, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>2){
		printf("sigmaP 1+2+3: ");
		for(int i=0;i<astringcount;i++){
			printf("%lf ",sigmaP[i]);
		}
		printf("\n");
	}
}

//compute sigma^Q for a given state up to a given number of alpha and beta terms in the trial wavefunction
void get_factored_sigmaQ(double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI,int print){
	int bstringcount = nchoosek(nmo,betae);
	get_factored_sigmaQBB(P,Q,Pprime,sigmaQ, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>1){
		printf("sigmaQ 2:     ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",sigmaQ[i]);
		}
		printf("\n");
	}

	get_factored_sigmaQAA(P,Q,Pprime,sigmaQ, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>0){
		printf("sigmaQ 1+2:   ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",sigmaQ[i]);
		}
		printf("\n");
	}
	get_factored_sigmaQAB(P,Q,Pprime,sigmaQ, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
	if(print>2){
		printf("sigmaQ 1+2+3: ");
		for(int i=0;i<bstringcount;i++){
			printf("%lf ",sigmaQ[i]);
		}
		printf("\n");
	}
}

//compute alpha-alpha component of sigma^P for given state
void get_factored_sigmaPAA(double *P,double *Q,double *Qprime, double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;
	if(debug_print){
		printf("PAA:\n");
	}
	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Iaprimestring = init_int_array(alphae);
	int *Iastring = init_int_array(alphae);
	int *Kastring = init_int_array(alphae);
	int *Jastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	double QB_prefactor = 0.0;
	//pre-compute sum_beta Q_beta*Q_beta
	QB_prefactor+=C_DDOT(bstringcount,Qprime,1,Q,1);
	if(fabs(QB_prefactor)>10E-6){	
		//loop up to highest excitation included in alpha wfn
		for(L=0;L<astringcount;L++){
			//loop over alpha 
			for(int kk=0;kk<alphae;kk++){
				Iastring[kk]=kk;
			}
			for(I=0;I<astringcount;I++){
				//first excitations
				for(l=0;l<alphae;l++){
					//one-electron coupling
					//only excite up to strings that can match Iaprime
					//for(k=0;k<nmo;k++){
					for(k=0;k<alphae;k++){
						//int sgnkl = excite(Iastring,k,Iastring[l],alphae,nmo,Kastring);
						int sgnkl = excite(Iastring,Iaprimestring[k],Iastring[l],alphae,nmo,Kastring);
						if(sgnkl != 0 ){
							int Kindex = stradr(Kastring,alphae,nmo);
							if(Kindex == L){
								//sigmaP[Kindex]+=sgnkl*mo_OEIprime[k][Iastring[l]]*P[I]*QB_prefactor;
								sigmaP[L]+=sgnkl*mo_OEIprime[Iaprimestring[k]][Iastring[l]]*P[I]*QB_prefactor;
								if(debug_print){
									printf("sigmaP[%d] += %d*h[%d][%d]*P[%d]*QB = %20.14lf\n",Kindex,sgnkl,Iaprimestring[k],Iastring[l],I,sigmaP[Kindex]);
								}
							}
						}
					}
					//two-electron coupling
					for(k=0;k<nmo;k++){
						int sgnkl = excite(Iastring,k,Iastring[l],alphae,nmo,Kastring);
						if(sgnkl != 0){
							//second excitations
							//for(i=0;i<betae;i++){
							for(i=0;i<alphae;i++){
								for(j=0;j<nmo;j++){
									//int sgnij = excite(Kastring,i,j,alphae,nmo,Jastring);
									int sgnij = excite(Kastring,Iaprimestring[i],j,alphae,nmo,Jastring);
									if(sgnij != 0 ){
										int Jindex = stradr(Jastring,alphae,nmo);
										if(Jindex == L){
											//sigmaP[Jindex]+=0.5*QB_prefactor*sgnkl*sgnij*mo_TEI[i*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Iastring[l]]*P[I];
											sigmaP[L]+=0.5*QB_prefactor*sgnkl*sgnij*mo_TEI[Iaprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Iastring[l]]*P[I];
											if(debug_print){
												printf("sigmaP[%d] += 1/2*%d*%d(%d%d|%d%d)*P[%d]*QB = %20.14lf\n",Jindex,sgnkl,sgnij,Iaprimestring[i],j,k,Iastring[l],I,sigmaP[Jindex]);
											}
										}
									}
								}
							}
						}
					}
				}
				next_combination(Iastring,nmo,alphae);	
			}
			next_combination(Iaprimestring,nmo,alphae);
		}
		if(debug_print){
			printf("\n");
		}
	}
	free(Iaprimestring);
	free(Iastring);
	free(Kastring);
	free(Jastring);
}

//compute beta-beta component of sigma^P
void get_factored_sigmaPBB(double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;
	if(debug_print){
		printf("PBB:\n");
	}

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibstring = init_int_array(betae);
	int *Ibprimestring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Jbstring = init_int_array(betae); 

	double Qfactor = 0.0;	

	for(i=0;i<betae;i++){
		Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	//loop over beta prime
	for(I=0;I<bstringcount;I++){
		//loop over beta
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		}
		for(J=0;J<bstringcount;J++){
			//first excitations
			for(l=0;l<betae;l++){
				//one-electron coupling
				for(k=0;k<betae;k++){
					int sgnkl = excite(Ibstring,Ibprimestring[k],Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						int Kindex = stradr(Kbstring,betae,nmo);
						if(Kindex == I){
						//	for(K=0;K<n_Aterms;K++){
						//	for(K=0;K<astringcount;K++){
							//	sigmaP[K]+=sgnkl*mo_OEIprime[Ibprimestring[k]][Ibstring[l]]*P[K]*Q[I]*Q[J];
							Qfactor+=sgnkl*mo_OEIprime[Ibprimestring[k]][Ibstring[l]]*Qprime[I]*Q[J];
						//	}

							//sigmaP[n_Aterms-1]+=sgnkl*mo_OEIprime[Ibprimestring[k]][Ibstring[l]]*P[n_Aterms-1]*Q[I]*Q[J];
							if(debug_print){
								//printf("sigmaP[%d] += %d*h[%d][%d]*P[%d]*Q[%d]*Q[%d] = %20.14lf\n",n_Aterms-1,sgnkl,Ibprimestring[k],Ibstring[l],n_Aterms-1,I,J,sigmaP[n_Aterms-1]);
							}
						}
					}
				}
				//two-electron coupling
				for(k=0;k<nmo;k++){
					int sgnkl = excite(Ibstring,k,Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						//second excitations
						for(i=0;i<betae;i++){
							for(j=0;j<nmo;j++){
								int sgnij = excite(Kbstring,Ibprimestring[i],j,betae,nmo,Jbstring);
								int Jindex = stradr(Jbstring,betae,nmo);
								if(Jindex == I && sgnij != 0){
									for(K=0;K<astringcount;K++){
										//sigmaP[K]+=0.5*sgnkl*sgnij*mo_TEI[Ibprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Ibstring[l]]*P[K]*Q[I]*Q[J];
									}
									Qfactor+=0.5*sgnkl*sgnij*mo_TEI[Ibprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Ibstring[l]]*Qprime[I]*Q[J];
									//sigmaP[n_Aterms-1]+=0.5*sgnkl*sgnij*mo_TEI[Ibprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Ibstring[l]]*P[n_Aterms-1]*Q[I]*Q[J];
									if(debug_print){
										//printf("sigmaP[%d] += 1/2*%d*%d*(%d%d|%d%d)*P[%d]*Q[%d]*Q[%d] = %20.14lf\n",n_Aterms-1,sgnkl,sgnij,Ibprimestring[i],j,k,Ibstring[l],n_Aterms-1,I,J,sigmaP[n_Aterms-1]);
									}
								}
							}
						}
					}
				}
			}

			next_combination(Ibstring,nmo,betae);	
		}
		next_combination(Ibprimestring,nmo,betae);	
	}
	for(i=0;i<astringcount;i++){
		sigmaP[i]+=P[i]*Qfactor;
	}
	if(debug_print){
		printf("\n");
	}
	free(Ibprimestring);
	free(Ibstring);
	free(Kbstring);
	free(Jbstring);
}

//compute alpha-beta component of sigma^P
void get_factored_sigmaPAB(double *P,double *Q,double *Qprime,double *sigmaP,int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;
	if(debug_print){
		printf("PAB:\n");
	}

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibprimestring = init_int_array(betae);
	int *Ibstring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Iastring = init_int_array(alphae); 
	int *Iaprimestring = init_int_array(alphae); 
	int *Kastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	for(i=0;i<betae;i++){
		Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	/*for(I=0;I<astringcount;I++){
		for(i=0;i<betae;i++){
			Ibprimestring[i]=i;
		}
		for(K=0;K<bstringcount;K++){
			for(i=0;i<betae;i++){
				Ibstring[i]=i;
			}
			for(L=0;L<bstringcount;L++){
				for(i=0;i<betae;i++){
					for(j=0;j<betae;j++){
						int sgnij = excite(Ibstring,Ibprimestring[i],Ibstring[j],betae,nmo,Kbstring);
						if(sgnij!=0){
						int Kbindex = stradr(Kbstring,betae,nmo);
						if(Kbindex == K){
							for(int kk=0;kk<alphae;kk++){
								Iastring[kk]=kk;
							}
							for(J=0;J<n_Aterms;J++){
								for(k=0;k<alphae;k++){
									for(l=0;l<alphae;l++){
										int sgnkl = excite(Iastring,Iaprimestring[k],Iastring[l],alphae,nmo,Kastring);
										if(sgnkl!=0){
										int Kaindex = stradr(Kastring,alphae,nmo);
										if(Kaindex == I){
											sigmaP[I]+=sgnkl*sgnij*Q[K]*Q[L]*P[I]*mo_TEI[Ibprimestring[i]*nmo*nmo*nmo+Ibstring[j]*nmo*nmo+Iaprimestring[k]*nmo+Iastring[l]];
											if(debug_print){
												printf("sigmaP[%d] += %d*%d*Q[%d]*Q[%d]*P[%d]*(%d%d|%d%d) = %20.14lf\n",I,sgnkl,sgnij,K,L,J,Ibprimestring[i],Ibstring[j],Iaprimestring[k],Iastring[l],sigmaP[I]); 
											}
										}}

									}
								}
								next_combination(Iastring,nmo,alphae);
							}
						}}
					}
				}
				next_combination(Ibstring,nmo,betae);
			}
			next_combination(Ibprimestring,nmo,betae);
		}
		next_combination(Iaprimestring,nmo,alphae);
	}*/
	//initialize beta ij density matrix
	double **beta_density = block_matrix(nmo,nmo);
	//loop over beta prime
	for(I=0;I<bstringcount;I++){	
		//loop over beta
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		} 
		for(J=0;J<bstringcount;J++){
			//excitation ij
			for(i=0;i<betae;i++){
				for(j=0;j<betae;j++){
					//only loop over i in the bra string and j in the ket string
					int sgnij = excite(Ibstring,Ibprimestring[i],Ibstring[j],betae,nmo,Kbstring);
					int Kindex = stradr(Kbstring,betae,nmo);
					if(Kindex == I && sgnij != 0){
						beta_density[Ibprimestring[i]][Ibstring[j]]+=sgnij*Qprime[I]*Q[J];
					}
				}
			}
			next_combination(Ibstring,nmo,betae);
		}
		next_combination(Ibprimestring,nmo,betae);
	}
	//print_mat(beta_density,nmo,nmo,outfile);
	//highest excitation for alpha'
	for(J=0;J<astringcount;J++){
		//loop over alpha
		for(int kk=0;kk<alphae;kk++){
			Iastring[kk]=kk;
		} 
		for(I=0;I<astringcount;I++){
			//excitations kl
			for(k=0;k<alphae;k++){
				for(l=0;l<alphae;l++){
					int sgnkl = excite(Iastring,Iaprimestring[k],Iastring[l],alphae,nmo,Kastring);
					if(sgnkl != 0 ){
						int Kindex = stradr(Kastring,alphae,nmo);
						if(Kindex == J ){
							for(i=0;i<nmo;i++){
								for(j=0;j<nmo;j++){
									sigmaP[J]+=beta_density[i][j]*sgnkl*mo_TEI[i*nmo*nmo*nmo+j*nmo*nmo+Iaprimestring[k]*nmo+Iastring[l]]*P[I];
									if(debug_print){
										printf("sigmaP[%d] += B[%d][%d]*%d*(%d%d|%d%d)*P[%d]\n",J,i,j,sgnkl,i,j,Iaprimestring[k],Iastring[l],I);
									}
								}
							}
						}
					}
				}
			}
			next_combination(Iastring,nmo,alphae);
		}
		next_combination(Iaprimestring,nmo,alphae);
	}
	if(debug_print){
		printf("\n");
	}
	free_block(beta_density);
	free(Ibprimestring);
	free(Ibstring);
	free(Kbstring);
	free(Iastring);
	free(Iaprimestring);
	free(Kastring);
}

//compute alpha-alpha component of sigma^Q
void get_factored_sigmaQAA(double *P,double *Q, double * Pprime, double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Iastring = init_int_array(alphae);
	int *Iaprimestring = init_int_array(alphae);
	int *Kastring = init_int_array(alphae);
	int *Jastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	//loop over alpha prime
	for(I=0;I<astringcount;I++){
		//loop over alpha
		for(int kk=0;kk<alphae;kk++){
			Iastring[kk]=kk;
		}
		for(J=0;J<astringcount;J++){
			//first excitations
			for(l=0;l<alphae;l++){
				//one-electron coupling
				for(k=0;k<alphae;k++){
					int sgnkl = excite(Iastring,Iaprimestring[k],Iastring[l],alphae,nmo,Kastring);
					if(sgnkl != 0 ){
						int Kindex = stradr(Kastring,alphae,nmo);
						if(Kindex == I){
							for(K=0;K<bstringcount;K++){
								sigmaQ[K]+=sgnkl*mo_OEIprime[Iaprimestring[k]][Iastring[l]]*Pprime[I]*P[J]*Q[K];
							}
							//sigmaQ[n_Bterms-1]+=sgnkl*mo_OEIprime[Iaprimestring[k]][Iastring[l]]*P[I]*P[J]*Q[n_Bterms-1];
							if(debug_print){
							//	printf("+ %d*%1.1f*h[%d][%d]Q[%d] ",sgnkl,Pa_term,Iaprimestring[k],Iastring[l],n_Bterms-1);
							}
						}
					}
				}
				//two-electron coupling
				for(k=0;k<nmo;k++){
					int sgnkl = excite(Iastring,k,Iastring[l],alphae,nmo,Kastring);
					if(sgnkl != 0){
						//second excitations
						for(i=0;i<alphae;i++){
							for(j=0;j<nmo;j++){
								int sgnij = excite(Kastring,Iaprimestring[i],j,alphae,nmo,Jastring);
								if(sgnij != 0){
									int Jindex = stradr(Jastring,alphae,nmo);
									if(Jindex == I){
										for(K=0;K<bstringcount;K++){
											sigmaQ[K]+=0.5*sgnkl*sgnij*mo_TEI[Iaprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Iastring[l]]*Pprime[I]*P[J]*Q[K];
										}
										//sigmaQ[n_Bterms-1]+=0.5*sgnkl*sgnij*mo_TEI[Iaprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Iastring[l]]*P[I]*P[J]*Q[n_Bterms-1];
										if(debug_print){
											//printf("+ 0.5*%d*%d*%1.1f*(%d %d|%d %d)Q[%d] ",sgnkl,sgnij,Pa_term,Iaprimestring[i],j,k,Iastring[l],n_Bterms-1);
										}
									}
								}
							}
						}
					}
				}
			}
			next_combination(Iastring,nmo,alphae);	
		}
		next_combination(Iaprimestring,nmo,alphae);	
	}
	if(debug_print){
		printf("\n");
	}
	free(Iastring);
	free(Iaprimestring);
	free(Kastring);
	free(Jastring);
}


//compute beta-beta component of sigma^Q
void get_factored_sigmaQBB(double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibprimestring = init_int_array(betae);
	int *Ibstring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Jbstring = init_int_array(betae); 
	for(i=0;i<betae;i++){
		//Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	double PA_prefactor = 0.0;
	//pre-compute sum_beta Q_beta*Q_beta
	PA_prefactor=C_DDOT(astringcount,Pprime,1,P,1);
	//PA_prefactor=1.0;
	//printf("PA prefactor = %lf\n",PA_prefactor);
	//loop up to highest excitation included in beta wfn
	for(L=0;L<bstringcount;L++){
		//loop over beta
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		}
		for(I=0;I<bstringcount;I++){
			//first excitations
			for(l=0;l<betae;l++){
				//one-electron coupling
				for(k=0;k<betae;k++){
					int sgnkl = excite(Ibstring,Ibprimestring[k],Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						int Kindex = stradr(Kbstring,betae,nmo);
						if(Kindex == L){
							sigmaQ[L]+=PA_prefactor*sgnkl*mo_OEIprime[Ibprimestring[k]][Ibstring[l]]*Q[I];
							/*if(debug_print && L == 0){
								printf("sigmaQ[%d] += %lf*%d*h[%d][%d]*Q[%d] ",L,PA_prefactor,sgnkl,Ibprimestring[k],Ibstring[l],I);
								for(int ii=0;ii<betae;ii++){
									printf("%d ",Ibprimestring[ii]);
								}
								printf("<%d' %d> ",Ibprimestring[k],Ibstring[l]);
								for(int ii=0;ii<betae;ii++){
									printf("%d ",Ibstring[ii]);
								}
								printf("\n");
							}*/
						}								
					}
				}
				//two-electron coupling
				for(k=0;k<nmo;k++){
					int sgnkl = excite(Ibstring,k,Ibstring[l],betae,nmo,Kbstring);
					if(sgnkl != 0){
						//second excitations
						for(i=0;i<betae;i++){
							for(j=0;j<nmo;j++){
								int sgnij = excite(Kbstring,Ibprimestring[i],j,betae,nmo,Jbstring);
								if(sgnij != 0){
									int Jindex = stradr(Jbstring,betae,nmo);
									if(Jindex == L){
										sigmaQ[L]+=PA_prefactor*0.5*sgnkl*sgnij*mo_TEI[Ibprimestring[i]*nmo*nmo*nmo+j*nmo*nmo+k*nmo+Ibstring[l]]*Q[I];
										/*if(debug_print && L == 0 && k == 5){
											printf("sigmaQ[%d] += 1/2*%lf*%d*%d*(%d%d|%d%d)*Q[%d] ",L,PA_prefactor,sgnkl,sgnij,Ibprimestring_ex[i],j,k,Ibstring_ex[l],I);
											for(int ii=0;ii<betae;ii++){
												printf("%d ",Ibprimestring[ii]);
											}
											printf("<%d' %d %d' %d> ",Ibprimestring[i],j,k,Ibstring[l]);
											for(int ii=0;ii<betae;ii++){
												printf("%d ",Ibstring[ii]);
											}
											printf("\n ");
											for(int ii=0;ii<betae;ii++){
												printf("%d ",Jbstring[ii]);
											}
											printf("[%d] <%d' %d > ",Jindex,Ibprimestring_ex[i],j);
											for(int ii=0;ii<betae;ii++){
												printf("%d ",Kbstring[ii]);
											}
											printf("[%d] <%d' %d > ",stradr(Kbstring,betae,nmo),k,Ibstring_ex[l]);
											for(int ii=0;ii<betae;ii++){
												printf("%d ",Ibstring[ii]);
											}
											printf("|%d %d ",nmo,betae);
											printf("\n");
										}*/
									}
								}
							}
						}
					}
				}
			}
			next_combination(Ibstring,nmo,betae);	
		}
		next_combination(Ibprimestring,nmo,betae);
	}
	if(debug_print){
		printf("\n");
	}
	free(Ibstring);
	free(Ibprimestring);
	free(Kbstring);
	free(Jbstring);
}

//compute alpha-beta component of sigma^Q
void get_factored_sigmaQAB(double *P,double *Q,double *Pprime,double *sigmaQ, int alphae, int betae, int nmo, double **mo_OEIprime, double * mo_TEI){

	int debug_print = 0;

	//dummy variables
	int i,j,k,l,I,J,K,L;
	//# strings
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae); 
	//initial strings
	int *Ibprimestring = init_int_array(betae);
	int *Ibstring = init_int_array(betae);
	int *Kbstring = init_int_array(betae);
	int *Iastring = init_int_array(alphae); 
	int *Iaprimestring = init_int_array(alphae); 
	int *Kastring = init_int_array(alphae); 
	for(i=0;i<alphae;i++){
		Iastring[i]=i;
		Iaprimestring[i]=i;
	}
	for(i=0;i<betae;i++){
		Ibstring[i]=i;
		Ibprimestring[i]=i;
	}
	//loop over alpha prime	
	//initialize alpha ij density matrix
	double **alpha_density = block_matrix(nmo,nmo);
	//loop over alpha prime
	for(I=0;I<astringcount;I++){	
		//loop over alpha 
		for(int kk=0;kk<alphae;kk++){
			Iastring[kk]=kk;
		}
		for(J=0;J<astringcount;J++){
			//excitation ij
			for(i=0;i<alphae;i++){
				for(j=0;j<alphae;j++){
					//only loop over i in the bra string and j in the ket string
					int sgnij = excite(Iastring,Iaprimestring[i],Iastring[j],alphae,nmo,Kastring);
					int Kindex = stradr(Kastring,alphae,nmo);
					if(Kindex == I && sgnij != 0){
						alpha_density[Iaprimestring[i]][Iastring[j]]+=sgnij*Pprime[I]*P[J];
					}
				}
			}
			next_combination(Iastring,nmo,alphae);
		}
		next_combination(Iaprimestring,nmo,alphae);
	}
	//highest excitation for beta'
	for(L=0;L<bstringcount;L++){
		//loop over beta 
		for(int kk=0;kk<betae;kk++){
			Ibstring[kk]=kk;
		}
		for(I=0;I<bstringcount;I++){
			//excitations kl
			for(k=0;k<betae;k++){
				for(l=0;l<betae;l++){
					int sgnkl = excite(Ibstring,Ibprimestring[k],Ibstring[l],betae,nmo,Kbstring);
					int Kindex = stradr(Kbstring,betae,nmo);
					if(Kindex == L && sgnkl != 0 ){
						for(i=0;i<nmo;i++){
							for(j=0;j<nmo;j++){
								sigmaQ[L]+=alpha_density[i][j]*sgnkl*mo_TEI[i*nmo*nmo*nmo+j*nmo*nmo+Ibprimestring[k]*nmo+Ibstring[l]]*Q[I];
								if(debug_print){
								//	printf("+ %1.1f*%d*(%d %d|%d %d)Q[%d] ",alpha_density[i][j],sgnkl,i,j,Ibprimestring[k],Ibstring[l],I);
								}
							}
						}
					}
				}
			}
			next_combination(Ibstring,nmo,betae);
		}
		next_combination(Ibprimestring,nmo,betae);
	}
	if(debug_print){
		printf("\n");
	}
	free_block(alpha_density);


	free(Ibstring);
	free(Ibprimestring);
	free(Kbstring);
	free(Iastring);
	free(Iaprimestring);
	free(Kastring);
}


	//get_factored_sigmaP(state, sigmaP, P, Q, n_Aterms, n_Bterms, alphae, betae, nmo, mo_OEIprime,  mo_TEI);

double get_energyP(int state, double **P, double **sigmaP, int n_Aterms, int alphae, int nmo){
	int astringcount = nchoosek(nmo,alphae);
	double E = C_DDOT(astringcount,&(sigmaP[state][0]),1,&(P[state][(n_Aterms-1)*astringcount]),1);
	return E;
}

double get_energyQ(int state, double **Q, double **sigmaQ, int n_Bterms, int betae, int nmo){
	int bstringcount = nchoosek(nmo,betae);
	double E = C_DDOT(bstringcount,&(sigmaQ[state][0]),1,&(Q[state][(n_Bterms-1)*bstringcount]),1);
	return E;
}

//dE/dP_alpha = 2*(sigma_alpha - E*P_alpha)
//make sure the newest addition to the wavefunction is normalized wrt previous iterations
void get_energy_gradientP(int state, double E, double **gradientP, double **sigmaP, double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI){
	int astringcount = nchoosek(nmo,alphae);
	E*=-1;
	C_DCOPY(astringcount,&(sigmaP[state][0]),1,&(gradientP[state][0]),1);
	C_DAXPY(astringcount,E,&(P[state][(n_Aterms-1)*astringcount]),1,&(gradientP[state][0]),1);
	C_DSCAL(astringcount,2,&(gradientP[state][0]),1);
}

//dE/dP_alpha = 2*(sigma_alpha - E*P_alpha)
//make sure the newest addition to the wavefunction is normalized wrt previous iterations
void get_energy_gradientQ(int state, double E,double **gradientQ, double **sigmaQ, double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEIprime, double *mo_TEI){
	int bstringcount = nchoosek(nmo,betae);
	E*=-1;
	C_DCOPY(bstringcount,&(sigmaQ[state][0]),1,&(gradientQ[state][0]),1);
	C_DAXPY(bstringcount,E,&(P[state][(n_Bterms-1)*bstringcount]),1,&(gradientQ[state][0]),1);
	C_DSCAL(bstringcount,2,&(gradientQ[state][0]),1);
}

void get_approx_diag_invHessP(int state, double *hessDiagP, double E, double **sigmaP, double **gradientP, int n_Aterms, int alphae, int nmo){

	int astringcount = nchoosek(nmo,alphae);
	for(int i=0;i<astringcount;i++){
		hessDiagP[i]=sigmaP[state][i]*gradientP[state][i]+E;
	}
	C_DSCAL(astringcount,-2,&(gradientP[state][0]),1);
}

//Davidson diagonalization for P table, holding Q fixed
int davidP(int state, int M, double *total_energy, double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI, double **mo_OEIprime, double *mo_TEI,int use_guess,int print)
{
	//dummy variables
	int i, j, k, L, I;
	double minimum;	
	int min_pos, numf, iter, *conv, converged, maxdim, skip_check;
	int *small2big, init_dim;
	int smart_guess =1;
	double *Adiag, **b, **bnew, **sigma, **G;
	double *lambda, **alpha, **f, *lambda_old;
	double norm, denom, diff;
	double BIGNUM = 10E100;

	//maximum number of iterations
	int MAXIT = 100;
	
	//convergence criterion for the ENERGY
	double cutoff = 10E-6;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);

	//dimension includes c0 normalization parameter
	int N = astringcount+1;

	//maximum Davidson subspace dimension (M is the number of roots sought)
	maxdim = 5 * M;

	b = block_matrix(maxdim, N);  /* current set of guess vectors,
				   stored by row */
	bnew = block_matrix(2*M, N); /* guess vectors formed from old vectors,
				stored by row*/
	sigma = block_matrix(maxdim,N); /* sigma vectors, stored by row*/
	G = block_matrix(maxdim, maxdim); /* Davidson mini-Hamitonian */
	f = block_matrix(maxdim, N); /* residual eigenvectors, stored by row */
	alpha = block_matrix(maxdim, maxdim); /* eigenvectors of G */
	lambda = init_array(maxdim); /* eigenvalues of G */
	lambda_old = init_array(maxdim); /* approximate roots from previous
				      iteration */
	Adiag = init_array(N);

	if(smart_guess) { /* Use eigenvectors of a sub-matrix as initial guesses */
		//set initial dimension
		if(N > maxdim){
			init_dim = (maxdim-1)*M;
		}
		else{ 
			init_dim = 2*M;
		}
		//if this witch is on, copy the previous best guess for current P
		if(use_guess){
			C_DCOPY(astringcount,&(P[state][(n_Aterms-1)*astringcount]),1,&(b[0][1]),1);	
		}
		//build sub-Hamiltonian element-by-element
		double **initG = block_matrix(init_dim,init_dim);
		for(i=0; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				initG[i][j] = build_single_Hamiltonian_element((n_Aterms-1)*astringcount+i,(n_Aterms-1)*astringcount+j,mo_OEI,mo_TEI,alphae,betae,nmo);

			}
		}
		fprintf(outfile,"init G P\n");
		print_mat(initG,init_dim,init_dim,outfile);

		//diagonalize sub-Hamiltonian
		sq_rsp(init_dim, init_dim, initG, lambda, 1, alpha, 1e-12);
		
		//set the first element of b (c0) to 1, then normalize
		int guess_index = 0;
		if(use_guess){
			guess_index = 1;
		}
		double c0 = 0.0;
		if(n_Aterms>1){
			c0 = 1.0;
			b[0][0] = c0;
			normalize(b[0],N);
		}
		for(i=guess_index; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				b[i][j+1] = alpha[j][i];
			}
			b[i][0]=c0;
			normalize(b[i],N);
		}
		
		//schmidt orthonormalize b vectors wrt converged P tables
		if(n_Aterms>1){
			for(i=0;i<n_Aterms-1;i++){
				for(j=0;j<init_dim;j++){
					double proj = C_DDOT(astringcount,&(b[j][1]),1,&(P[state][i*astringcount]),1);
					C_DAXPY(astringcount,-proj,&(P[state][i*astringcount]),1,&(b[j][1]),1);
				}
				normalize(&(b[i][1]),astringcount);
			}
		}
		
		//normalize first b vector
		normalize(b[0],N);
		
		//schmidt orthonormalize b vectors with respect to each other
		for(i=1;i<init_dim;i++){
			for(j=0;j<i;j++){
				double proj = C_DDOT(N,b[j],1,b[i],1);
				C_DAXPY(N,-proj,b[j],1,b[i],1);
			}
			normalize(b[i],N);
		}
		
		free_block(initG);
	}

	//set current dimension L to initial dimension
	L = init_dim;
	iter =0;
	converged = 0;
	conv = init_int_array(M); /* boolean array for convergence of each
			       root */
	//keep track of current number of sigma vectors in memory
	int built = 0;

	//build diagonal elements of H
	for(i=0;i<N;i++){
		Adiag[i] = build_single_Hamiltonian_element(i*astringcount+n_Bterms-1,i*astringcount+n_Bterms-1,mo_OEI,mo_TEI,alphae,betae,nmo);
	}
	//ITERATE
	while(converged < M && iter < MAXIT) {
		//keep track of whether this iteration required subspace collapse
		skip_check = 0;
		if(print){
			 printf("iter = %d\n", iter); 
		}
		//print guess vectors (b are guesses)
		fprintf(outfile,"B %d (L = %d)\n",iter,L);
		print_mat(b,L,N,outfile);
	
		//check orthonormality of b vectors
		int print_overlap = 1;
		if(print_overlap){
			double **overlap = block_matrix(L,L);
			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(b[0][0]), N, 0.0, &(overlap[0][0]), L);
			fprintf(outfile,"P overlap\n");
			print_mat(overlap,L,L,outfile);
			free_block(overlap);	
		}
		//build sigma vectors
		if(n_Aterms>1){
			for(i=built;i<L;i++){
				get_augmented_sigmaP(&(b[i][0]), &(Q[state][(n_Bterms-1)*bstringcount]), P, Q, &(sigma[i][0]), total_energy[0], n_Aterms, n_Bterms, alphae, betae,  nmo, mo_OEIprime, mo_TEI, 1);
				built++;
			}
		}
		else{
			for(i=built;i<L;i++){
				printf("b %d: ",iter);
				for(j =0;j<astringcount;j++){
					printf("%lf ",b[i][j+1]);
				}
				printf("\n");
				get_factored_sigmaP(&(b[i][1]),&(Q[state][0]),&(Q[state][0]),&(sigma[i][1]), alphae, betae, nmo, mo_OEIprime,  mo_TEI,0);
				built++;
			}
		}

		//print sigma vectors
		fprintf(outfile,"sigma P %d\n",iter);
		print_mat(sigma,built,N,outfile);
	
		/* form mini-matrix */
		C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(sigma[0][0]), N, 0.0, &(G[0][0]), maxdim);

		//print quadratic forms
		fprintf(outfile,"subHP %d\n",iter);
		print_mat(G,L,L,outfile);

		/* diagonalize mini-matrix */
		sq_rsp(L, L, G, lambda, 1, alpha, 1e-12);
		fprintf(outfile,"P eigenvalues\n");
		for(j=0;j<M;j++){
			fprintf(outfile,"%lf\n",lambda[j]);
		}
		fprintf(outfile,"subHP eigenvectors\n");
		print_mat(alpha,L,L,outfile);

		/* form preconditioned residue vectors */
		for(k=0; k < M; k++) {//rows
			for(I=0; I < N; I++) { //cols
				f[k][I] = 0.0;
				for(i=0; i < L; i++) {
					f[k][I] += alpha[i][k] * (sigma[i][I] - lambda[k] * b[i][I]);
				}
				denom = lambda[k] - Adiag[I];
				if(fabs(denom) > 10e-6) {
					f[k][I] /= denom;
				}
				else{
					f[k][I] = 0.0;
				}
			}
		}

		/* normalize each residual */
		for(k=0; k < M; k++) {
			norm = 0.0;
			for(I=0; I < N; I++) {
				norm += f[k][I] * f[k][I];
			}
			norm = sqrt(norm);
			for(I=0; I < N; I++) {
				if(norm > 1e-6) {
					f[k][I] /= norm;
				}
				else {
					f[k][I] = 0.0;
				}
			}
		}

		/* schmidt orthogonalize the f[k] against the set of b[i] and add
		new vectors */
		for(k=0,numf=0; k < M; k++){
			if(schmidt_add(b, L, N, f[k])) { 
				printf("added\n");L++; numf++; 
			}
		}

		/* If L is close to maxdim, collapse to two guesses per root */
		if(maxdim - L < M) {
			if(print) {
				printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
				printf("Collapsing eigenvectors.\n");
				fprintf(outfile,"P COLLAPSE ==============\n");
			}
			for(i=0;i<L;i++){
				for(j=0;j<N;j++){
					sigma[i][j]=0.0;
				}
			}
			for(i=0; i < 2*M; i++) {
				memset((void *) bnew[i], 0, N*sizeof(double));
				for(j=0; j < L; j++) {
					for(k=0; k < N; k++) {
						bnew[i][k] += alpha[j][i] * b[j][k];
						//bnew[i][k] += G[j][i] * b[j][k];
					}
				}
			}
			/* copy new vectors into place */
			for(i=0; i < 2*M; i++){ 
				for(k=0; k < N; k++){
					b[i][k] = bnew[i][k];
				}
			}
			if(n_Aterms==1){
				for(i=0;i<M;i++){
					b[i][0]=0.0;
					normalize(&(b[i][0]),N);
				}
			}
			skip_check = 1;
			built = 0;
			L = 2*M;
		}

		/* check convergence on all roots */
		if(!skip_check) {
			converged = 0;
			zero_int_array(conv, M);
			if(print) {
				printf("Root      Eigenvalue       Delta  Converged?\n");
				printf("---- -------------------- ------- ----------\n");
			}
			for(k=0; k < M; k++) {
				diff = fabs(lambda[k] - lambda_old[k]);
				if(diff < cutoff) {
					conv[k] = 1;
					converged++;
				}
				lambda_old[k] = lambda[k];
				if(print) {
					printf("%3d  %20.14f %4.3e    %1s\n", k, lambda[k], diff,
					 conv[k] == 1 ? "Y" : "N");
				}
			}
		}

		iter++;
	}

	/* generate final eigenvalues and eigenvectors */
	if(converged == M) {
	//copy converged energies
	for(i=0;i<M;i++){
		total_energy[i] = lambda[i];
	}
	double **v = block_matrix(N,M);
		for(i=0; i < M; i++) {
		//eps[i] = lambda[i];
			for(j=0; j < L; j++) {
				for(I=0; I < N; I++) {
					v[I][i] += alpha[j][i] * b[j][I];
				}
			}
		}
	//printf("\nc0 P = %lf\n\n",v[0][state]);
		if(print) printf("Davidson algorithm converged in %d iterations.\n", iter);
		//scale previously converged vectors by c0
		if(n_Aterms>1){
			C_DSCAL((n_Aterms-1)*astringcount,v[0][state],&(P[state][0]),1);
		}
		//copy newest converged vector
		for(I=0;I<astringcount;I++){
			P[state][(n_Aterms-1)*astringcount+I]=v[I+1][state];
		}
		//normalize newest converged vector
		normalize(&(P[state][(n_Aterms-1)*astringcount]),astringcount);
		printf("Converged energy = %lf\n",total_energy[0]);
	fprintf(outfile,"||||||||||Converged eigenvectors P|||||||||||||\n");
	print_mat(v,N,M,outfile);
	fprintf(outfile,"|||||||||||||||||||||||||||||||||||||||||||||||\n");
	}

  free(conv);
  free_block(b);
  free_block(bnew);
  free_block(sigma);
  free_block(G);
  free_block(f);
  free_block(alpha);
  free(lambda);
  free(lambda_old);

  return converged;
}


//same as davidP but it forms sigma^Q vectord
//see davidP for more extensive comments
int davidQ(int state, int M,double *total_energy, double **P, double **Q, int n_Aterms, int n_Bterms, int alphae, int betae, int nmo, double **mo_OEI,double **mo_OEIprime, double *mo_TEI,int use_guess,int print){
	int i, j, k, L, I;
	double minimum;
	int min_pos, numf, iter, *conv, converged, maxdim, skip_check;
	int *small2big, init_dim;
	int smart_guess =1;
	double *Adiag, **b, **bnew, **sigma, **G;
	double *lambda, **alpha, **f, *lambda_old;
	double norm, denom, diff;
	double BIGNUM = 10E100;
	int MAXIT = 100;
	double cutoff = 10E-6;
	int astringcount = nchoosek(nmo,alphae);
	int bstringcount = nchoosek(nmo,betae);
	int N = bstringcount+1;
	maxdim = 5 * M;

	b = block_matrix(maxdim, N);  /* current set of guess vectors,
				   stored by row */
	bnew = block_matrix(2*M, N); /* guess vectors formed from old vectors,
				stored by row*/
	sigma = block_matrix(maxdim,N); /* sigma vectors, stored by row*/
	G = block_matrix(maxdim, maxdim); /* Davidson mini-Hamitonian */
	f = block_matrix(maxdim, N); /* residual eigenvectors, stored by row */
	alpha = block_matrix(maxdim, maxdim); /* eigenvectors of G */
	lambda = init_array(maxdim); /* eigenvalues of G */
	lambda_old = init_array(maxdim); /* approximate roots from previous
				      iteration */
		Adiag = init_array(N);
	if(n_Bterms == 1){
		smart_guess = 1;
		printf("Q Sub-Hamiltonian Guess\n");
	}
	if(smart_guess) { /* Use eigenvectors of a sub-matrix as initial guesses */
		if(N > maxdim){
			init_dim = (maxdim-1)*M;
		}
		else{ 
			init_dim = 2*M;
		}
		if(use_guess){
			C_DCOPY(bstringcount,&(Q[state][(n_Bterms-1)*bstringcount]),1,&(b[0][1]),1);	
		}
		
		//build sub-Hamiltonian element-by-element
		double **initG = block_matrix(init_dim,init_dim);
			for(i=0; i < init_dim; i++) {
				for(j=0; j < init_dim; j++){
					initG[i][j] = build_single_Hamiltonian_element(i*astringcount+n_Bterms-1,j*astringcount+n_Bterms-1,mo_OEI,mo_TEI,alphae,betae,nmo);
				}
			}
		
		fprintf(outfile,"init G Q\n");
		print_mat(initG,init_dim,init_dim,outfile);
		
		//diagonalize sub-Hamiltonian
		sq_rsp(init_dim, init_dim, initG, lambda, 1, alpha, 1e-12);
		
		//set c0 = 1 and normalize b vectors
		int guess_index=0;
		if(use_guess){
			guess_index = 1;
		}
		double c0 = 0.0;
		if(n_Aterms>1){
			c0 = 1.0;
			b[0][0] = c0;
			normalize(b[0],N);
		}
		for(i=guess_index; i < init_dim; i++) {
			for(j=0; j < init_dim; j++){
				b[i][j+1] = alpha[j][i];
			}
			b[i][0]=c0;
			normalize(b[i],N);
		}
		//schmidt ON b wrt converged Q tables
		if(n_Bterms>1){
			for(i=0;i<n_Bterms-1;i++){
				for(j=0;j<init_dim;j++){
					double proj = C_DDOT(bstringcount,&(b[j][1]),1,&(Q[state][i*bstringcount]),1);
					C_DAXPY(bstringcount,-proj,&(Q[state][i*bstringcount]),1,&(b[j][1]),1);
				}
				normalize(&(b[i][1]),bstringcount);
			}
		}
		//normalize first b vector
		normalize(b[0],N);
		//schmidt ON b vectors wrt each other
		for(i=1;i<init_dim;i++){
			for(j=0;j<i;j++){
				double proj = C_DDOT(N,b[j],1,b[i],1);
				C_DAXPY(N,-proj,b[j],1,b[i],1);
			}
			normalize(b[i],N);
		}
		

		free_block(initG);
	}

	L = init_dim;
	iter =0;
	converged = 0;
	conv = init_int_array(M); /* boolean array for convergence of each
			       root */
	int built = 0;

	//diagonal elements of H
	for(i=0;i<bstringcount;i++){
		Adiag[i] = build_single_Hamiltonian_element(astringcount*(n_Aterms-1)+i,astringcount*(n_Aterms-1)+i,mo_OEI,mo_TEI,alphae,betae,nmo);
	}

	//ITERATE
	while(converged < M && iter < MAXIT) {

		skip_check = 0;
		if(print){
			 printf("\niter = %d\n", iter); 
		}
		fprintf(outfile,"B Q %d\n",iter);
		print_mat(b,L,N,outfile);
		//check orthonormality of b
		int print_overlap = 1;
		if(print_overlap){
			double **overlap = block_matrix(L,L);
			C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(b[0][0]), N, 0.0, &(overlap[0][0]), L);
			fprintf(outfile,"Q overlap\n");
			print_mat(overlap,L,L,outfile);
			free_block(overlap);	
		}
		//build sigma vectors
		if(n_Bterms>1){
			for(i=built;i<L;i++){
				get_augmented_sigmaQ(&(P[state][(n_Aterms-1)*astringcount]), &(b[i][0]), P, Q, &(sigma[i][0]), total_energy[0], n_Aterms, n_Bterms, alphae, betae,  nmo, mo_OEIprime, mo_TEI, 1);
				built++;
			}
		}
		else{
			for(i=built;i<L;i++){
				printf("b %d: ",iter);
				for(int kk =1;kk<bstringcount+1;kk++){
					printf("%lf ",b[i][kk]);
				}
				printf("\n");
				double *sigref = &(sigma[i][1]);
				get_factored_sigmaQ(&(P[state][0]),&(b[i][1]),&(P[state][0]),sigref,  alphae, betae, nmo, mo_OEIprime,  mo_TEI,0);
				built++;
			}
		}
		fprintf(outfile,"SIGMA Q %d\n",iter);
		print_mat(sigma,built,N,outfile);
		/* form mini-matrix */
		C_DGEMM('n','t', L, L, N, 1.0, &(b[0][0]), N,&(sigma[0][0]), N, 0.0, &(G[0][0]), maxdim);
		fprintf(outfile,"subHQ %d\n",iter);
		print_mat(G,L,L,outfile);
		/* diagonalize mini-matrix */
		sq_rsp(L, L, G, lambda, 1, alpha, 1e-12);
		fprintf(outfile,"Q eigenvalues\n");
		for(int kk=0;kk<M;kk++){
			fprintf(outfile,"%lf\n",lambda[kk]);
		}
		fprintf(outfile,"subHQ eigenvectors\n");
		print_mat(alpha,L,L,outfile);
		/* form preconditioned residue vectors */
		for(k=0; k < M; k++) {//rows
			for(I=0; I < N; I++) { //cols
				f[k][I] = 0.0;
				for(i=0; i < L; i++) {
					f[k][I] += alpha[i][k] * (sigma[i][I] - lambda[k] * b[i][I]);
				}
				denom = lambda[k] - Adiag[I];
				if(fabs(denom) > 1e-6) {
					f[k][I] /= denom;
				}
				else{
					f[k][I] = 0.0;
				}
			}
		}

		/* normalize each residual */
		for(k=0; k < M; k++) {
			norm = 0.0;
			for(I=0; I < N; I++) {
				norm += f[k][I] * f[k][I];
			}
			norm = sqrt(norm);
			for(I=0; I < N; I++) {
				if(norm > 1e-6) {
					f[k][I] /= norm;
				}
				else {
					f[k][I] = 0.0;
				}
			}
		}

		/* schmidt orthogonalize the f[k] against the set of b[i] and add
		new vectors */
		for(k=0,numf=0; k < M; k++){
			if(schmidt_add(b, L, N, f[k])) { 
				printf("added\n");L++; numf++; 
			}
		}

		/* If L is close to maxdim, collapse to one guess per root */
		if(maxdim - L < M) {
			if(print) {
				printf("Subspace too large: maxdim = %d, L = %d\n", maxdim, L);
				printf("Collapsing eigenvectors.\n");
				fprintf(outfile,"Q COLLAPSE ==============\n");
			}
			for(i=0;i<L;i++){
				for(j=0;j<N;j++){
					sigma[i][j]=0.0;
				}
			}
			for(i=0; i < 2*M; i++) {
				memset((void *) bnew[i], 0, N*sizeof(double));
				for(j=0; j < L; j++) {
					for(k=0; k < N; k++) {
						bnew[i][k] += alpha[j][i] * b[j][k];
						//bnew[i][k] += G[j][i] * b[j][k];
					}
				}
			}
			/* copy new vectors into place */
			for(i=0; i < 2*M; i++){ 
				for(k=0; k < N; k++){
					b[i][k] = bnew[i][k];
				}
				//normalize(&(b[i][0]),N);
			}
			if(n_Bterms==1){
				for(i=0;i<M;i++){
					b[i][0]=0.0;
					normalize(&(b[i][0]),N);
				}
			}
			skip_check = 1;
			built = 0;
			L = 2*M;
		}

		/* check convergence on all roots */
		if(!skip_check) {
			converged = 0;
			zero_int_array(conv, M);
			if(print) {
				printf("Root      Eigenvalue       Delta  Converged?\n");
				printf("---- -------------------- ------- ----------\n");
			}
			for(k=0; k < M; k++) {
				diff = fabs(lambda[k] - lambda_old[k]);
				if(diff < cutoff) {
					conv[k] = 1;
					converged++;
				}
				lambda_old[k] = lambda[k];
				if(print) {
					printf("%3d  %20.14f %4.3e    %1s\n", k, lambda[k], diff,
					 conv[k] == 1 ? "Y" : "N");
				}
			}
		}

		iter++;
	}
	/* generate final eigenvalues and eigenvectors */
	if(converged == M) {
	//copy covnerged energies
	for(i=0;i<M;i++){
		total_energy[i] = lambda[i];
	}
	double **v = block_matrix(N,M);
		for(i=0; i < M; i++) {
		//eps[i] = lambda[i];
			for(j=0; j < L; j++) {
				for(I=0; I < N; I++) {
					v[I][i] += alpha[j][i] * b[j][I];
				}
			}
		}
	//printf("\nc0 Q = %lf\n\n",v[0][state]);
		if(print) printf("Davidson algorithm converged in %d iterations.\n", iter);
		//scale previously converged tables by c0 
		if(n_Bterms>1){
			C_DSCAL((n_Bterms-1)*bstringcount,v[0][state],&(Q[state][0]),1);
		}
		//copy and normalize newest converged tables
		for(I=0;I<bstringcount;I++){
			Q[state][(n_Bterms-1)*bstringcount+I]=v[I+1][state];
			Q[state][I]=v[I+1][state];
		}
		normalize(&(Q[state][(n_Bterms-1)*bstringcount]),bstringcount);
	fprintf(outfile,"||||||||||||||||||||||Eigenvectors Q||||||||||||||||||||\n");
	print_mat(v,N,M,outfile);
	fprintf(outfile,"||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	}

  free(conv);
  free_block(b);
  free_block(bnew);
  free_block(sigma);
  free_block(G);
  free_block(f);
  free_block(alpha);
  free(lambda);
  free(lambda_old);

  return converged;
}

}}


	/*Adiag = init_array(N);
		for(i=0; i < N; i++) { 
			Adiag[i] = build_single_Hamiltonian_element(i*astringcount+n_Bterms-1,i*astringcount+n_Bterms-1,mo_OEIprime,mo_TEI,alphae,betae,nmo);
			printf("diag(%d) <%d|H|%d> = %lf\n",i,i*astringcount+n_Bterms-1,i*astringcount+n_Bterms-1,Adiag[i]);
		}
		for(i=0; i < init_dim; i++) {
			minimum = Adiag[0];
			min_pos = 0;
			for(j=1; j < N; j++){
				if(Adiag[j] < minimum) { 
					minimum = Adiag[j]; min_pos = j; 
				}
			}
			b[i][min_bpos] = 1.0; 
			Adiag[min_pos] = BIGNUM; 
			lambda_old[i] = minimum;
		}*/



	/*double *DIAG = init_array(N);
	for(int i=0;i<astringcount;i++){
		DIAG[i]=1.0;//sqrt(astringcount);
	}
	int *IPRINT = init_int_array(2);
	double EPS = 10e-3;
	double XTOL = 10e-16;
	double *W = init_array(N*(2*M+1)+2*M);
	int IFLAG = 0; 
	printf("LBFGS\n");
	int iter = 0;
		double *G = init_array(N);
	while(iter<2){
		printf("iter %d\n",iter);
		IPRINT[0] = 1;
		IPRINT[1] = 3;
		for(int i=0;i<astringcount;i++){
			sigmaP[state][i]=0.0;
		}
		for(int i=0;i<bstringcount;i++){
			sigmaQ[state][i]=0.0;
		}
		//get_factored_sigmaP(state, sigmaP, P, Q, 1, 1, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
		get_factored_sigmaQ(state, sigmaQ, P, Q, 1, 1, alphae, betae, nmo, mo_OEIprime,  mo_TEI);
		double F = get_energyQ(0,Q,sigmaQ, 1, betae, nmo);	
		double f = get_energyP(0,P,sigmaP, 1, alphae, nmo);	
		printf("E = %lf\n",F);
		printf("E = %lf\n",f);
		//double *G = init_array(N);
		get_energy_gradientP(0,F, &G, sigmaP, P, Q, 1, 1, alphae, betae, nmo, mo_OEIprime, mo_TEI);
		double *hessDiagP = init_array(N);
		get_approx_diag_invHessP(state, hessDiagP, F, sigmaP, &G, 1, alphae, nmo);
		lbfgs_(&N, &M, &(P[state][0]), &f, G, &DIAGCO, hessDiagP, IPRINT, &EPS, &XTOL, W, &IFLAG);
		free(hessDiagP);
		for(int i=0;i<astringcount;i++){
			printf("%lf\n",P[state][i]);
		}
		//normalize(&(P[state][0]),astringcount);
		printf("IFLAG = %d\n",IFLAG);
		if(IFLAG==1){IFLAG = 0;}
		//free(G);
		iter++;
	}*/
	/*free_block(sigmaP);
	free_block(sigmaQ);

	free_block(P);
	free_block(Q);*/


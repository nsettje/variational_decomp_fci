#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#ifndef _psi_plugins_variational_decomp_fci_MOint_h
#define _psi_plugins_variational_decomp_fci_MOint_h

using namespace boost;

namespace psi { namespace variational_decomp_fci {

void initialize_MO_constants(int *nmo, int *alphae, int *betae, double *nuc_rep_energy,double *eSCF, const char * molname, const char * basisname);
void MO_transform(double **mo_OEI, double *mo_TEI, int nmo);
void OEItrans(double **S, double **C, double **mo_OEI,int nmo);
void TEItrans(double *mo_TEI, double **C);

void initialize_MO_constants(int *nmo, int *alphae, int *betae,double *nuc_rep_energy, double *eSCF, const char * molname, const char *basisname){
	// Molecule and wavefunction objects
	shared_ptr<Molecule> mol = Process::environment.molecule();                                                              
	std::string molname_ = mol->name();
	molname = molname_.c_str();
	shared_ptr<Wavefunction> wfn = Process::environment.reference_wavefunction();
	// Basis set object		
	shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
	shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, mol, "BASIS");
																  
	// Integral factory
	shared_ptr<IntegralFactory> integral = shared_ptr<IntegralFactory>(new IntegralFactory(aoBasis,aoBasis,aoBasis,aoBasis));
	std::string basisname_ = aoBasis->name();
	basisname = basisname_.c_str();
	*nmo = wfn->nmo();
	*alphae = wfn->nalpha();
	*betae = wfn->nbeta();
	*nuc_rep_energy = mol->nuclear_repulsion_energy();
	*eSCF = wfn->reference_energy();
}
void MO_transform(double **mo_OEI, double *mo_TEI, int nmo){
	// Molecule and wavefunction objects
	shared_ptr<Molecule> mol = Process::environment.molecule();                                                              
	shared_ptr<Wavefunction> wfn = Process::environment.reference_wavefunction();
	// Basis set object		
	shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
	shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, mol, "BASIS");

	shared_ptr<IntegralFactory> integral = shared_ptr<IntegralFactory>(new IntegralFactory(aoBasis,aoBasis,aoBasis,aoBasis));
	// Integral
	shared_ptr<TwoBodyAOInt> eri = shared_ptr<TwoBodyAOInt>(integral->eri(0));
	//matrix factory
	int* nsopi = wfn->nsopi();
	shared_ptr<MatrixFactory> factory = shared_ptr<MatrixFactory>(new MatrixFactory);
	factory->init_with(1, nsopi, nsopi);
																  
	// MO coefficients
	SharedMatrix _C = wfn->Ca();
	double **C = _C->pointer();
	//printf("Built HF C matrix\n"); 
																  
	//AO overlap
	shared_ptr<OneBodyAOInt> Sint(integral->ao_overlap());
	shared_ptr<Matrix> sMat(factory->create_matrix("Overlap"));
	//Sint->compute(sMat);
	//S = sMat->pointer();
	shared_ptr<OneBodyAOInt> tOBI(integral->ao_kinetic());
	shared_ptr<OneBodyAOInt> vOBI(integral->ao_potential());
	//printf("Built S matrix\n");
	shared_ptr<Matrix> tMat(factory->create_matrix("Kinetic"));
	shared_ptr<Matrix> vMat(factory->create_matrix("Potential"));
	shared_ptr<Matrix> HCore(factory->create_matrix("T+V"));
																  
			//printf("here\n");													  
	vOBI->compute(vMat);
			//printf("there\n");													  
	tOBI->compute(tMat);
																  
	HCore->copy(tMat);
	HCore->add(vMat);
	double **ao_OEI = HCore->pointer();
	//printf("Built core Hamiltonian\n");
																  
	//Important constants of the system
	OEItrans(mo_OEI, C, ao_OEI, nmo);      
	TEItrans(mo_TEI,C);                
}


//Transform OEI to MO basis (C*OEI*C)
void OEItrans(double **mo_OEI, double **C, double **ao_OEI, int nmo){                      
	double **temp = block_matrix(nmo,nmo);
	C_DGEMM('t','n',nmo,nmo,nmo,1.0,C[0],nmo,ao_OEI[0],nmo,0.0,temp[0],nmo); //temp=C*S
	C_DGEMM('n','n',nmo,nmo,nmo,1.0,temp[0],nmo,C[0],nmo,0.0,mo_OEI[0],nmo); //mo_OEI=temp*C
	print_mat(mo_OEI,nmo,nmo,outfile);
	free_block(temp);
}

//Computes TEI and transforms them to MO basis
void TEItrans(double *mo_TEI, double **C){
	//molecule, wavefunction, basis, and integral objects                                                      
	//[could probably move these but makes this call easier]
	shared_ptr<Molecule> mol = Process::environment.molecule();
	shared_ptr<Wavefunction> wfn =
			Process::environment.reference_wavefunction();	
	shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
	shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, mol, "BASIS");
	int nmo = wfn->nmo();
	int nshell=aoBasis->nshell();
	// Integral factory
	shared_ptr<IntegralFactory> integral = shared_ptr<IntegralFactory>(new
			IntegralFactory(aoBasis,aoBasis,aoBasis,aoBasis));
	// Integral
	shared_ptr<TwoBodyAOInt> eri = shared_ptr<TwoBodyAOInt>(integral->eri(0));
	
	double *ao_TEI=init_array(nmo*nmo*nmo*nmo);
				  
	const double *buffer=eri->buffer();
	for(int MU=0, uvpoI=0; MU<aoBasis->nshell(); MU++) {
    int nummu = aoBasis->shell(MU).nfunction();
    for(int NU=0; NU<aoBasis->nshell(); NU++) {
      int numnu = aoBasis->shell(NU).nfunction();
      for(int RHO=0; RHO<aoBasis->nshell(); RHO++) {
        int numrho = aoBasis->shell(RHO).nfunction();
        for(int SIG=0; SIG<aoBasis->nshell(); SIG++) {
          int numsig = aoBasis->shell(SIG).nfunction();
          // Compute a shell quartet
          eri->compute_shell(MU,NU,RHO,SIG);
          // Loop over all the integrals in the shell quartet.
          for(int mu=0, munu=0, index=0; mu<nummu; mu++) {
            int omu = aoBasis->shell(MU).function_index() + mu;
            for(int nu=0; nu<numnu; nu++, munu++) {
              int onu = aoBasis->shell(NU).function_index() + nu;
              for(int rho=0; rho<numrho; rho++) {
                int orho = aoBasis->shell(RHO).function_index() + rho;
                for(int sig=0; sig<numsig; sig++, index++,uvpoI++) {
                  int osig = aoBasis->shell(SIG).function_index() + sig;
                  ao_TEI[(((omu*nmo+onu)*nmo+orho)*nmo)+osig]=buffer[index];
                  //printf("<%d %d | %d %d> = %f\n",onu,omu,orho,osig,buffer[index]);
                }
              }
            }
          }
        } // End SIGMA Shell
      } // End RHO Shell
    } // End NU Shell
  } // End MU Shell

	//Transform TEI to MO basis 
	// MO basis TEI stored in mo_TEI (uvpo->pqrs)
	// SIGMA -> S
	double *tempS=init_array(nmo*nmo*nmo*nmo);
	for(int mu=0;mu<nmo;mu++){
	  for(int nu=0;nu<nmo;nu++){
	    for(int rho=0;rho<nmo;rho++){
	      for(int sig=0;sig<nmo;sig++){
	      	for(int s=0;s<nmo;s++){
	        	tempS[(((mu*nmo+nu)*nmo+rho)*nmo)+s]+=ao_TEI[(((mu*nmo+nu)*nmo+rho)*nmo)+sig]*C[sig][s];
	        }
	      }
	    }
	  }
	}
	free(ao_TEI);
	// RHO -> R
	double *tempR=init_array(nmo*nmo*nmo*nmo);
	for(int mu=0;mu<nmo;mu++){
	  for(int nu=0;nu<nmo;nu++){
	  	for(int rho=0;rho<nmo;rho++){
	  		for(int s=0;s<nmo;s++){
	        	for(int r=0;r<nmo;r++){
	        	tempR[(((mu*nmo+nu)*nmo+r)*nmo)+s]+=tempS[(((mu*nmo+nu)*nmo+rho)*nmo)+s]*C[rho][r];
	        }
	      }
	    }
	  }
	}
	free(tempS);
	// NU -> Q
	double *tempQ=init_array(nmo*nmo*nmo*nmo);
	for(int mu=0;mu<nmo;mu++){
		for(int nu=0;nu<nmo;nu++){
			for(int r=0;r<nmo;r++){
				for(int s=0;s<nmo;s++){
					for(int q=0;q<nmo;q++){
	        	tempQ[(((mu*nmo+q)*nmo+r)*nmo)+s]+=tempR[(((mu*nmo+nu)*nmo+r)*nmo)+s]*C[nu][q];
	        }
	      }
	    }
	  }
	}
	free(tempR);
	// MU -> P
	for(int mu=0;mu<nmo;mu++){
		for(int q=0;q<nmo;q++){
			for(int r=0;r<nmo;r++){
			  for(int s=0;s<nmo;s++){
				for(int p=0;p<nmo;p++){
	        	mo_TEI[(((p*nmo+q)*nmo+r)*nmo)+s]+=tempQ[(((mu*nmo+q)*nmo+r)*nmo)+s]*C[mu][p];
	        }
	      }
	    }
	  }
	}
	free(tempQ);
	for(int p=0;p<nmo;p++){
	for(int q=0;q<nmo;q++){
	for(int r=0;r<nmo;r++){
	for(int s=0;s<nmo;s++){
		//fprintf(outfile,"(%d %d | %d %d) = %20.14lf\n",p,q,r,s,mo_TEI[p*nmo*nmo*nmo+q*nmo*nmo+r*nmo+s]);
	}
	}
	}
	}
	
}

}}

#endif

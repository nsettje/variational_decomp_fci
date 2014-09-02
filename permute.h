#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#ifndef _psi_plugins_variational_decomp_fci_permute_h
#define _psi_plugins_variational_decomp_fci_permute_h

using namespace boost;
namespace psi { namespace variational_decomp_fci {

int nchoosek(int n, int k);
void Hklprime(double *mo_TEI, double **mo_OEI, double **hklprime,int nmo);
int excite(int *Istring, int k, int l, int nelec, int nmo, int *excited);     
int excite_fast(int *Istring, int k, int l, int nelec, int nmo, int *excited);
int zindex(int k, int l, int nelec, int nmo);
int stradr(int *str, int nelec, int nmo);
int trunc_stradr(int *str, int nelec, int nmo);
unsigned int next_combination(int * const ar, int n, unsigned int k);
unsigned int next_combination_char(unsigned char * ar, int n,unsigned int k);
int next_combination_fast(int * const ar, int n, int k);

//! Binomial coefficient n!/(k!*(n-k)!) [probably slow]
int nchoosek(int n, int k){
	//check that we are not choosing more objects than exist
        if(k>n){
            return 0;
        }
        
        int nfact=1;
	for(int i=n-k+1;i<n+1;i++){
		nfact*=i;
	}
	int kfact=1;
	for(int i=1;i<k+1;i++){
		kfact*=i;
	}
	
	int denom = nfact/kfact;
	return denom;
}

int excite(int *Istring, int k, int l, int nelec, int nmo, int *excited){
	int i,j;
	int sign;
	int *sort=init_int_array(nelec);
	//C_DCOPY(nelec,Istring,1,excited,1);
	if(k==l){
		for(i=0;i<nelec;i++){
			if(Istring[i]==l){
				for(j=0;j<nelec;j++){
					excited[j] = Istring[j];
				}
				sign = 1;
				return sign;
			}
		}
	}
	int ann_ind = -1;
	for(i=0;i<nelec;i++){
		if(Istring[i]==k){
			sign = 0;
			return sign;
		}
		if(Istring[i]==l){
			ann_ind = i;
		}
		sort[i] = Istring[i];
	}
	if(ann_ind<0){
		return 0;
	}
	//printf("ann = %d |",ann_ind);
	sort[ann_ind] = k;
	int min,min_pos;
	for(i=0;i<nelec;i++){
		min = sort[0];
		for(j=0;j<nelec;j++){
			if(sort[j]<=min){
				min = sort[j];
				excited[i] = sort[j];
				min_pos = j;
				//sign_count += j - i;
			}
		}
		sort[min_pos] = INT_MAX;
	}
	int sign_count = 0;
	for(i=0;i<nelec;i++){
		if(Istring[i]<k && Istring[i] != l){
			sign_count++;
		}
	}
	sign_count+=ann_ind;
	//printf("sign_count = %d |",sign_count);
	if(sign_count % 2 == 0){
		sign = 1;
		return sign;
	}
	else{
		sign = -1;
		return sign;
	}
	free(sort);
	/*for(i=0;i<nelec;i++){
		if(Istring[i] == l){
			if(k == l){
				sign = 1;
				return sign;
			}
			excited[i] = Istring[i];
		}	
	}*/
	//if k = l, copy ground to excited
	/*if(k==l){
		
	}
	//check if ground string contains the annihilation index l
	int match_index=-1;
	for(i=0;i<nelec;i++){
		if(Istring[i]==l){
			match_index = i;
			Istring[i] = -1;
		}
	}
	if(match_index<0){
		return 0;	
	}
	//check if ground string contains the creation index k
	int create_index;
	for(i=0;i<nelec;i++){
		if(Istring[i]==k){
			Istring[match_index] = l;
			return 0;
		}
		if(Istring[i]<k && Istring[i]>0){
			create_index = i;
		}
	}
	Istring[match_index]=l;
	for(i=0;i<create_index+1;i++){
		excited[i] = Istring[i];
	}
	excited[create_index+1] = k;
	for(i=create_index+2;i<*/
	
	/*int sign;
	//for excitations into same level, skip most of the work by checking if
	//the string contains the excitation index
	if(k==l){
		for(int j=0;j<nelec;j++){
			if(l==Istring[j]){
				for(int i=0;i<nelec;i++){
					excited[i]=Istring[i];
				}
				sign=1;
				return(sign);
			}
		}
	}
	
	int ann_index=0;

	//loop until annihilation matches an index
	for(int i=0;i<nelec;i++){
		if(l==Istring[i]){
			break;
		}
		ann_index++;
	}
	//if no index matches, sign = 0
	if(ann_index==nelec){
		return(0);
	}
	
	int create_sweep=0;
	//loop until creation matches an index
	for(int i=0;i<nelec;i++){
		if(k==Istring[i]){
			break;
		}
		create_sweep++;
	}
	//if creation matches an index, sign = 0
	if(create_sweep!=nelec){
		return(0);
	}
	
	for(int i=0;i<nelec;i++){
		excited[i]=Istring[i];
	}
		excited[ann_index]=k;
	int perm=ann_index;
	int sign_perm=ann_index;
	int buff;
	if(k<l){
		while(excited[ann_index]<excited[ann_index-1 && ann_index>0]){
			buff=excited[ann_index];
			excited[ann_index]=excited[ann_index-1];
			excited[ann_index-1]=buff;
			ann_index--;
			perm++;
		}
	}
	else if(k>l){
		while(excited[ann_index]>excited[ann_index+1] && ann_index+1<nelec){
			buff=excited[ann_index];
			excited[ann_index]=excited[ann_index+1];
			excited[ann_index+1]=buff;
			ann_index++;
			perm++;
		}
	}
	sign_perm+=perm;
	if(sign_perm%2==0){
		sign=1;
	}
	else{
		sign=-1;
	}
		
	return sign;*/	
	
}

//given a combination of length k over n possible elements,
//compute the next combination in lexical ordering
unsigned int next_combination_char(unsigned char * ar, int n,unsigned int k){
    unsigned int finished = 0;
    unsigned int changed = 0;
    unsigned int i;
    unsigned int place;

    for (i = k - 1, place = 1; !finished && !changed; i--,place++) {
        if (ar[i] < n - place) {
            //Increment this element 
            ar[i]++;
            if (i < k - 1) {
                /* Make the elements after it the same */
                int j;
                for (j = i + 1; j < k ; j++) {
                    ar[j] = ar[j - 1]+1;
                }
            }
            changed = 1;
        }
        finished = i == 0;
    }
    return changed;
}


//given a combination of length k over n possible elements,
//compute the next combination in lexical ordering
unsigned int next_combination(int * const ar, int n,unsigned int k){
    unsigned int finished = 0;
    unsigned int changed = 0;
    unsigned int i;
    unsigned int place;

    for (i = k - 1, place = 1; !finished && !changed; i--,place++) {
        if (ar[i] < n - place) {
            //Increment this element 
            ar[i]++;
            if (i < k - 1) {
                /* Make the elements after it the same */
                int j;
                for (j = i + 1; j < k ; j++) {
                    ar[j] = ar[j - 1]+1;
                }
            }
            changed = 1;
        }
        finished = i == 0;
    }
    return changed;
}

// finds K-H z index [equation 11 from Knowles, Handy. Chem. Phys. Lett. Vol. 111 1984 pgs. 315-321] 
int zindex(int k, int l, int nelec, int nmo){
	int z=0;
	//valence
//	printf("z_%d = ",k);
	if(k==nelec-1){
		z=l-nelec+1;
	//	printf("%d ",z);
		}
	else if(l>k){
	//core
		for(int m=nmo-l;m<nmo-k;m++){
			z+=nchoosek(m,nelec-k-1)-nchoosek(m-1,nelec-k-2);
		//	printf("+ (%d,%d) - (%d,%d) ",m,nelec-k-1,m-1,nelec-k-2);
		}
	//	printf("= %d",z);
	}
//	printf("\n");
	return z;
}

//finds string address from sum over z indices. Make sure the orbitals actually exist...
//[equation 12 from Knowles, Handy. Chem. Phys. Lett. Vol. 111 1984 pgs. 315-321] 
int stradr(int *str, int nelec,int nmo){
	int adr=0;
	for(int i=0;i<nelec;i++){
		adr+=zindex(i,str[i],nelec,nmo);
	}
	return adr;
}

//finds the lexical ordering index of a string that lives in the space spanned by all excitations up to the highest level minus one
//use this for truncated CI
//use stradr for FCI
int trunc_stradr(int *str, int nelec, int nmo){
	int adr=0;
	int factor = nmo-nelec;
	for(int i = 0;i<nelec-1;i++){
		adr+=factor*(str[i]-i);
	}
	adr+=str[nelec-1]-nelec+1;
	return adr;
}

}}
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#endif

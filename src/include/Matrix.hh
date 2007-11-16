
#ifndef _MATRIX_HH
#define _MATRIX_HH

#include "Vector.hh"

#include<iostream>

namespace replaceleda{

    //matrix class
    template <class T>
    class mmatrix{
    public:
	mmatrix() : N(0), M(0){};

	mmatrix(std::vector< mvector<T> > vv){
	    //this one is not ment...
	    mx.clear();
	    for(uint i = 0; i < vv.size(); ++i)
		mx.push_back(vv[i]);
	    N=vv.size();
	    M=vv[0].size();
	
	    //use this one
	/*    mx.clear();
	    M=vv.size();
	    N=vv[0].size();
	    for(uint m=0 ; m < M ; ++m){
		mvector<T> dummy(N);
		for(uint n=0; n<N; ++n)
		    dummy[n] = vv[m][n];
		mx.push_back(dummy);
	    }
	*/
	}
	
	mmatrix(uint n, uint m) : N(n), M(m) {
	    mvector<T> D(n*m);
	    uint c=0;
	    for(uint i=0;i<n;++i){
		mvector<T> dummy;
		for(uint j=0;j<m;++j)
		    dummy.push_back(D[c++]);
		mx.push_back(dummy);
	    }
	}
	
	mmatrix(uint n, uint m, mvector<T> D) : N(n), M(m){
	    if(n*m <= D.size()){
		uint c=0;
		for(uint i=0;i<n;++i){
		    mvector<T> dummy;
		    for(uint j=0;j<m;++j)
			dummy.push_back(D[c++]);
		    mx.push_back(dummy);
		}
	    }
	}
	
	virtual ~mmatrix(){ 
	    for(uint i=0;i<N;++i)
		mx[i].clear();
	    mx.clear();
	}

	T elem(uint n, uint m){ return mx[n][m]; }

	void set(uint n, uint m, T elem){
	    mx[n].set(m, elem);
	}
	    
	//overloaad std operations
	mvector<T> operator[](uint i) const { return mx[i];}
	mvector<T>& operator[](uint i) {return mx[i];}

	mvector<T> row(uint i) const { return mx[i]; }
	mvector<T>& row(uint i){ return mx[i]; }

	mvector<T> col(uint j){
	    mvector<T> result;
	    for(int i=0;i<dim1();++i)
		result.push_back(mx[i][j]);
	    return result;
	}

	void operator=(mmatrix<T> m){
	    mx.clear();
	    for(int j=0;j<m.dim1();++j)
		mx.push_back(m[j]);
	    N=m.dim1();
	    M=m.dim2();
	}

	bool operator==(mmatrix<T> m){
	    if(m.dim1() != dim1() or m.dim2() != dim2())
		return false;
	    for(uint i=0;i<dim1();++i)
		for(uint j=0;j<dim2();++j)
		    if(m[i][j] != mx[i][j])
			return false;
	    return true;
	}

	bool operator!=(mmatrix<T> m){
	    return !operator==(m);
	}
	
	mmatrix<T> operator*(T r){
	    mvector<T> dummy;	
	    for(uint i=0;i<N;++i)
		for(uint j=0;j<M;++j)
		    dummy.push_back(r*mx[i][j]);
	    mmatrix<T> result(N,M,dummy);
	    return result;	
	}

	mmatrix<T> operator/(T r){
	    return operator*( T(1/r));
	}

	mmatrix<T> operator+(mmatrix<T> m){
	    mvector<T> dummy;
	    if(m.dim1() == dim1() && m.dim2() == dim2()){
		for(uint i=0;i<N;++i)
		    for(uint j=0;j<M;j++)
			dummy.push_back(m[i][j] + mx[i][j]);
		mmatrix<T> result(N,M,dummy);
		return result;
	    }
	    mmatrix<T> result(0,0);
	    return result;
	}

	mmatrix<T> operator-(mmatrix<T> m){
	    return operator+(m*T(-1));
	}
	
	void operator*=(T r){
	    for(uint i=0;i<N;++i)
		for(uint j=0;j<M;++j)
		    mx[i][j]*=r;
	}
	void operator/=(T r){ operator*=(T(1/r));}
	
	//matrix operations
	int dim1(){ return N; }
	int dim2(){ return M; }
	
	T operator() (uint i, uint j) const { return mx[i][j];}

	T& operator() (uint i, uint j){ return mx[i][j];}
	
	mmatrix<T> trans(){
	    std::vector<mvector<T> > dummy;
	    for(int j=0;j<dim2();++j)
		dummy.push_back(col(j));
	    
	    mmatrix<T> result(dummy);
	    return result;
	}
	
	void print(std::ostream &os){
	    for(int i=0;i<dim1();++i){
		for(int j=0;j<dim2();++j)
		    os << mx[i][j] << " ";
		os << std::endl;	    
	    }
	}
	void read(std::istream &is){
	    for(uint i=0;i<dim1();++i){
		for(uint j=0;j<dim2();++j)
		    is >> mx[i][j];
	    }
	}
    protected:
	mvector < mvector<T> > mx;
	uint N,M;
	
    };
/*
    template <class T>
    std::ostream& operator<< (std::ostream &os, mmatrix<T> m){
	for(uint i=0;i<m.dim1();++i){
	    for(uint j=0;j<m.dim2();++j)
		//os << m.elem(i,j) << " ";
		os << m[i][j] << " ";
	    os << std::endl;	    
	}
	return os;
    }

    template <class T>
    std::istream& operator>> (std::istream &is, mmatrix<T> &m){
	uint N,M;
	T e;
	mvector<T> dummy;
	std::string line;
	
	//read N and M
	is >> N; is >> M;
	while(is >> e)
	    dummy.push_back(e);
	if(dummy.dim() == N * M){
	    mmatrix<T> result(N,M,dummy);
	    m = result;	
	}
	return is;
    }
    
    template <class T>
    mmatrix<T> transpose(mmatrix<T> m){
	return m.trans();
    }
*/    
/*
    template <class T>
    std::ostream& operator<< (std::ostream &os, mmatrix<T> m);

    template <class T>
    std::istream& operator>> (std::istream &is, mmatrix<T> &m);

    template <class T>
    mmatrix<T> transpose(mmatrix<T> m);
*/
/*
    std::ostream& operator<< (std::ostream &os, mmatrix<int> m);
    std::ostream& operator<< (std::ostream &os, mmatrix<double> m);

    std::istream& operator>> (std::istream &is, mmatrix<int> &m);
    std::istream& operator>> (std::istream &is, mmatrix<double> &m);

    mmatrix<int> transpose(mmatrix<int> m);
    mmatrix<double> transpose(mmatrix<double> m);
*/
}

#endif //_Matrix_hh

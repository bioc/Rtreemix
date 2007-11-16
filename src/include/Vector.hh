#ifndef _VECTOR_HH
#define _VECTOR_HH

#include <iostream>
#include <vector>
#include <cmath>

namespace replaceleda{

    template <class T>
    class mvector{
    public:
	mvector(): D(0) { v.resize(0);}
	mvector(uint _dim): D(_dim){ v.resize(_dim); }
	mvector(T a, T b): D(2){
	    v.push_back(a); v.push_back(b);
	}
	mvector(T a, T b, T c): D(3){
	    v.push_back(a); v.push_back(b); v.push_back(c);
	}
	virtual ~mvector(){v.clear();}

	//vector standard operations
	T operator[](uint i) const { return v[i];}

	T& operator[](uint i) {return v[i];}

	void set(uint i, T val){ v[i]=val;}

	// calculations
	void operator=(mvector<T> v2){
	    v.clear();
	    D = v2.dim();
	    for(uint i=0;i<D;++i)
		v.push_back(v2[i]);	    
	}

	mvector<T> operator*(T r){
	    mvector<T> result;
	    for(uint i=0;i<D;++i)
		result.push_back(r*v[i]);
	    return result;
	}

	mvector<T> operator/(T r){
	    return operator*(T(1/r));
	}

	mvector<T> operator+(mvector<T> v2){
	    mvector<T> result;
	    for(uint i=0;i<D;++i)
		result.push_back(v[i] + v2[i]);
	    return result;
	}

	mvector<T> operator-(mvector<T> v2){
	    mvector<T> result;
	    for(uint i=0;i<D;++i)
		result.push_back(v[i]-v2[i]);
	    return result;

	}

	void operator*=(T r){
	    for(uint i=0;i<D;++i)
		v[i] *= r;
	}
	void operator/=(T r){
	    for(uint i=0;i<D;++i)
		v[i] /= r;
	}
	void operator+=(mvector<T> v2){
	    for(uint i=0;i<D;++i)
		v[i] +=v2[i];
	}
	void operator-=(mvector<T> v2){
	    for(uint i=0;i<D;++i)
		v[i] -=v2[i];
	}
	

	//comparision
	bool operator==(mvector<T> v2){
	    if (dim() != v2.dim())
		return false;
	    for(uint i=0;i<D;++i)
		if(v[i] != v2[i])
		    return false;	    
	    return true;
	}

	bool operator!=(mvector<T> v2){
	    return ! operator==(v2);
	}
	
	//std::vector operations
	void push_back(T elem){
	    v.push_back(elem);
	    D = v.size();
	}
	void clear(){
	    v.clear();
	    D = 0;
	}
	uint size(){
	    return D;
	}

	//vector specific operations
	int dim(){ return D;}
	
	T operator*(mvector<T> v2){
	    T result = T(0);
	    for(uint i=0;i<D;++i)
		result += v[i]*v2[i];
	    return result;
	}
	
	double sqr_length(){
	    T result = 0.0;
	    for(uint i=0;i<D;++i)
		result += v[i]*v[i];
	    return double(result);
	}

	double angle(mvector<T> v2){
	    T num = operator*(v2);
	    double denom = length() * v2.length();
	    double degree = acos ((double) num / denom);
	    if (degree > 180.0 )
		return 360 - degree;
	}
	
	double length(){
	    return sqrt(sqr_length());	    
	}

    protected:
	uint D;
	std::vector<T> v;
    };

/*
    template <class T>
    std::ostream& operator<< (std::ostream &os, mvector<T> v);

    template <class T>
    std::istream& operator>> (std::istream &is, mvector<T> &v);
*/
/*
    template <class T>
    std::ostream& operator<< (std::ostream &os, mvector<T> v){
	for(uint j=0;j<v.size();++j)
	    os << v[j] << " ";
	os << std::endl;
	return os;
    } 

    template <class T>
    std::istream& operator>> (std::istream &is, mvector<T> &v){
	T e;
	v.clear();
	while(is >> e)
	    v.push_back(e);	
	return is;
    }
*/
}


#endif // _Vector_hh

#ifndef _Array_hh
#define _Array_hh

#include <stdlib.h>
#include <vector>

namespace replaceleda{

template <class T>
class array : public std::vector<T> {
protected:
    typedef std::vector<T> backpoint;

public:
    array(uint num = 0) : backpoint(num) {}
    virtual ~array(){}
    
    int size() { return backpoint::size();}

    int low(){ return 0; }
    int high() { return (size() -1) ;}

    void permute(){
	std::vector<T> dummy;
	while(backpoint::size() > 0){
	    int r = int((backpoint::size() * (double) rand()) / (double) (RAND_MAX + 1.0));
	    dummy.push_back(backpoint::at(r));
	    backpoint::erase(backpoint::begin()+r);
	}
	backpoint::resize(dummy.size());
	for(uint i = 0; i < dummy.size(); ++i)
	    backpoint::at(i) = dummy[i];
    }

protected:

};


}


#endif //_Array_hh

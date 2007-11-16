#ifndef _SET_HH
#define _SET_HH

#include <set>

namespace replaceleda {

    template <class T>
    class set : public std::set<T> {
    protected:
	typedef std::set<T> backpoint;
    public:
	set () {};	
	virtual ~set() {}

	bool member(T item){
	    return (backpoint::find(item) != backpoint::end());
	}
	int size(){return backpoint::size();}
    };


}


#endif

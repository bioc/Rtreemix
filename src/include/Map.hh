#ifndef _MAP_HH
#define _MAP_HH

#include "List.hh"
#include <iterator>
#include <map>

namespace replaceleda{

    template <class T1, class T2>
    class map : public std::map<T1, T2> {
    protected:
	typedef std::map<T1, T2> backpoint;

    public:
	map() {}
	virtual ~map() {}

	bool defined(T1 index){ return backpoint::find(index) != backpoint::end();}
	list<T1> getDefined(){
	    list<T1> result;
	    typename backpoint::iterator it;
	    for(it = backpoint::begin() ; it != backpoint::end(); ++it)
	    	result.append(it->first);
	    return result;
	}
	
    protected:
	
    };

}

#endif

#ifndef _LIST_HH
#define _LIST_HH

#include <deque>
#include <list>
#include <vector>

namespace replaceleda {

    template <class T>
    class list : public std::deque<T> {	
    private:
	typedef std::deque<T> backpoint;
    public:
	list( uint num = 0, T item = T(0)) : backpoint(num, item) {}
	virtual ~list() {}
	
	void remove(T item){
	    uint i = 0;
	    while (i < backpoint::size() && backpoint::at(i) != item)
		++i;

	    if (i < backpoint::size())
		backpoint::erase( backpoint::begin() + i);	    
	}

	void append(T item){ backpoint::push_back(item); }
	T contents(uint pos) { return backpoint::at(pos); }

	uint max(int(&cmp)(const T&,const T&)){
	    uint mpos = 0;
	    
	    for(uint i = 0; i < backpoint::size(); ++i)
		if (cmp(backpoint::at(i), backpoint::at(mpos)) > 0)
		    mpos = i;    
	    return mpos;
	}
	uint min(int(&cmp)(const T&, const T&)){
	    uint mpos = 0;
	    for(uint i = 0; i < backpoint::size(); ++i)
		if (cmp(backpoint::at(i), backpoint::at(mpos)) < 0)
		    mpos = i;
	    return mpos;
	}

	uint length() { return backpoint::size(); }	

	T head(){ return backpoint::front();}

	T push(T item){ backpoint::push_front(item); return item; }
	T pop() { T result = backpoint::front(); backpoint::pop_front(); return result;}

	T get_item(int pos) {if ((uint)pos < backpoint::size())return backpoint::at(pos); return T(0);}
    
	void permute(){
	    std::vector<T> dummy;
	    while(backpoint::size() > 0){
		int r = int((backpoint::size() * (double) std::rand()) / (double) (RAND_MAX + 1.0));
		dummy.push_back(backpoint::at(r));
		backpoint::erase(backpoint::begin()+r);
	    }
	    backpoint::resize(dummy.size());
	    for(uint i = 0; i < dummy.size(); ++i)
		backpoint::at(i) = dummy[i];
	}

	void sort() {
	    std::list<T> dummy;
	    for (uint r = 0; r < length(); ++r)
		dummy.push_back(backpoint::at(r));
	    backpoint::clear();
	    dummy.sort();
	    while (dummy.size() > 0){
		backpoint::push_back(dummy.front());
		dummy.pop_front();
	    }
	}

	uint find(T item){
	    uint i = 0;
	    while( i < backpoint::size() and backpoint::at(i) != item)
		++i;

	    return i;	    
	}
        
	T next(T item, uint begin){
	    if (begin)
		return backpoint::at(0);
	    else {
		uint pos = find(item);
		if ((pos+1) < backpoint::size())
		    return backpoint::at(pos+1);
	    }
	    return item;
	}

	bool last(T item){
	    return item == backpoint::back();
	}

    private:	
    };
}



#endif //_LIST_HH

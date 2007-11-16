#ifndef _QUEUE_HH
#define _QUEUE_HH

#include <queue>

namespace replaceleda {

    template <class T>
    class queue : public std::queue<T> {
    protected:
	typedef std::queue<T> backpoint;
    public:
	queue(){}
	virtual ~queue(){}

	void append(T item){ backpoint::push(item);}
	T pop() { T result = backpoint::front(); backpoint::pop(); return result; }
    protected:
    };

    template<class T1, class T2>
    class pq_elem {
    public:
	T1 first;
	T2 second;
    public:
	pq_elem(){}
	pq_elem(T1 _first, T2 _second): first(_first), second(_second) {}
	virtual ~pq_elem(){};

	bool operator<(const pq_elem<T1, T2> &el) const {
	    return first > el.first; 
	}
    };

    typedef uint pq_item;

    template<class T1, class T2>
    class p_queue : public std::priority_queue<pq_elem<T1, T2> > {
    public:
	//typedef pq_elem<T1, T2> pq_item;
    protected:
	typedef std::priority_queue<pq_elem<T1, T2> > backpoint;
    public:
	p_queue(){}
	virtual ~p_queue(){}
        
	void insert(T1 opart, T2 spart){ pq_elem<T1, T2> nelem(opart, spart); backpoint::push(nelem);}
	//pq_elem<T1, T2> find_min(){ uint s = backpoint::size(); if(s>0) --s; return backpoint::at(s);}
	pq_item find_min() { return 0;}

	T2 inf(pq_item pos){ return backpoint::top().second;}

	void del_item(pq_item pos){ backpoint::pop();}
	
    };
}

#endif

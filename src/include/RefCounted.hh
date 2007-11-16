#ifndef _REFCOUNTED_HH
#define _REFCOUNTED_HH


#include <iostream>

namespace replaceleda { 

    class RefCounted {
    public:
	virtual ~RefCounted() {
	    //ASSERT(0 == m_nRefCount);
	    if (0 != m_nRefCount)
		std::cerr << "WARINING: # of counter to object is not zero!" << std::endl;
	}
	int getRefCount() { return m_nRefCount;}
    protected:
	RefCounted() {
	    m_nRefCount = 0;
	}
	// Never copy m_nRefCount
	RefCounted(RefCounted const& another) {
	    m_nRefCount = 0;
	}
	RefCounted& operator=(RefCounted const& another) {
	    return *this;
	}
	
	class ProtectConstructor { 
	public: 
	    ProtectConstructor(RefCounted const* pTarget) { 
		//ASSERT(m_pTarget); 
		m_pTarget->AddRef(); 
	    } 
	    ~ProtectConstructor() { 
		--m_pTarget->m_nRefCount; 
	    }
	private:
	    RefCounted const* m_pTarget;
	};
	
    private:
	mutable int m_nRefCount;
	
	void AddRef() const {
	    m_nRefCount += 1;
	}
	void DecRef() const {
	    m_nRefCount -= 1;
	}
	void Release() const {
	    DecRef();
	    if (0 == m_nRefCount) {
		//delete this;
		Destroy();
	    }
	}
	void Destroy() const {
	    delete this;
	}
	friend class RefCountPtrBase;
	friend class RefCounted::ProtectConstructor;
	
    };

    /**
     * Reference counting pointer
     **/
    class RefCountPtrBase 
    {
    public:
	RefCountPtrBase(RefCounted const* pTarget = 0): m_pTarget(pTarget) {
	    if (m_pTarget)
		m_pTarget->AddRef();
	}   
	
	RefCountPtrBase(RefCountPtrBase const& another): m_pTarget(another.m_pTarget) {
	    if (m_pTarget)
		m_pTarget->AddRef();
	}
	
	~RefCountPtrBase() {
	    if (m_pTarget)
		m_pTarget->Release();
	}      
	
	RefCountPtrBase& operator= (RefCounted const* pTarget) {
	    if (pTarget)
		pTarget->AddRef();
	    if (m_pTarget)
		m_pTarget->Release();
	    m_pTarget = pTarget;
	    return *this;        
	}
	
	RefCountPtrBase& operator= (RefCountPtrBase const& another) {
	    if (another.m_pTarget)
		another.m_pTarget->AddRef();
	    if (m_pTarget)
		m_pTarget->Release();
	    m_pTarget = another.m_pTarget;
	    return *this;        
	}
	
	RefCounted const& operator*() const {
	    return *m_pTarget;
	}
	
	RefCounted const* operator->() const {
	    return m_pTarget;
	}
    private:
	RefCounted const* m_pTarget;
    };


    /**
     * Typed reference counting pointers
     **/
    template<typename T>
    class RefCountPtr
    {
    public:
	RefCountPtr(T* pTarget = 0): m_pTarget(pTarget) {};
	RefCountPtr<T>& operator= (T * pTarget) {
	    m_pTarget = pTarget;
	    return *this;
	}
	T& operator*() const {
	    return (T&)(*m_pTarget);
	}
	T* operator->() const {
	    return (T*)(m_pTarget.operator->());
	}
	operator T* () const {
	    return operator->();
	}
    private:
	RefCountPtrBase m_pTarget;
    };
} // namespace


#endif //_REFCOUNTED_HH

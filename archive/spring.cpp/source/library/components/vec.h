// vector
template < typename T >
struct Vec : public vector<T>
{
	// constructor
	Vec() {}
	
	// fill
	Vec(int n)
	{
		vector<T>::resize(n);
	}
	
    // fill
	inline void resize(int n)
	{
        vector<T>::clear();
		vector<T>::resize(n);
	}
        	
	// fill
	inline void fill(T elem)
	{
		for (int i = 0; i < vector<T>::size(); i++)
			vector<T>::operator[](i) = elem;
	}
	
	inline T* ptr()
	{
        return &(*this)[0];
    }

	T operator[] (int i) const
	{
		assert (i < vector<T>::size());		
		return vector<T>::operator[](i);
	}
	
	T& operator[] (int i)
	{		
		//assert (i < vector<T>::size());
		return vector<T>::operator[](i);
	}
	
    inline void append(T element)
    {
        vector<T>::push_back(element);
    }
    
    inline int size()
    {
        return int(vector<T>::size());
    }    
};

// vector for bool
template < >
struct Vec < bool > : public vector<bool>
{
	// fill
	inline void fill(bool elem)
	{
		for (unsigned int i = 0; i < vector<bool>::size(); i++)
			vector<bool>::operator[](i) = elem;
	}
};

#include <iostream>
#include <map>

// hash
template < typename THash, typename TKey >
class Hash
{
	map < THash, TKey > list;
public:
	
	// constructor
	Hash()
	{
		list.clear();
	}
	
	// add key
	void insert (THash a, TKey b)
	{
		list.insert(make_pair (a, b));
	}
	
	// get key
	TKey get (THash a)
	{
           // get key
	    typename map < THash, TKey >::iterator mi = list.find(a);
	    if (mi == list.end())
               return -1;
        
           // return
           return mi->second;
	}	

	// get key
	int size ()
	{
           return list.size();
	}
	
	// get key
	void clear ()
	{
           list.clear();
	}	
};



#ifndef _ListMin
#define _ListMin

// list of minimal elements
template < typename TType >
struct ListMin
{
	// scores and list
	Vec < TType > 	litem;
	Vec < double > 	lscore;

	// size
	int nsize;

	// index of minimal value
	int index;

	// construct
	void construct(int nsize)
	{
        if (nsize < 1)
            Msg::error ("ListMin::construct()", "ListMin requires a size larger than zero.");

		// copy
		this->nsize = nsize;

		// items
		litem.resize(nsize);

		// scores
		lscore.resize(nsize);

		// initialize
		initialize();
	}

	// reset
	void initialize ()
	{
		lscore.fill(LARGE);
		index = 0;
	}

	// update
	inline void update (double score, TType item)
	{
		if (score < lscore[index])
		{
			// set new item
			lscore[index]	= score;
			litem[index]	= item;

			// search new maximal value and backup index
			for (int i = 0; i < nsize; i++)
				if (lscore[i] > lscore[index])
					index = i;
		}
	}

	// size
	int size()
	{
		return nsize;
	}

	// sort
	void sort()
	{
		Lib::sort(lscore, litem);
	}
};

#endif

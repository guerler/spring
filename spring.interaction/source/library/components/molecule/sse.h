// secondary structure element
struct SSE
{    
		int		start, end;
		int		n;
		char	type;
		
		// copy constructor
		SSE(SSE *sse)
        {
			this->start  = sse->start;
			this->end    = sse->end;
			this->type   = sse->type;
			this->n      = sse->n;
		}
		
		SSE(char type, int start = 0, int end = 0)
        {
			this->start  = start;
			this->end    = end;
			this->type   = type;
			n            = end - start + 1;
		}
		
        // recieve		
		inline int length() const 	{ return n; }
		inline int getStart() const { return start; }
        inline int getEnd() const 	{ return end; }
		inline char getType() const { return type;  }

        // alter		
		void setStart(int start)
        {
			this->start  = start;
			this->n      = end - start + 1;
		}
		
		void setEnd(int end)
        {
			this->end    = end;
			this->n      = end - start + 1;
		}			
};

	

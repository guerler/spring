// clock
class Clock
{
	// parameters
	time_t tmS, dur;	
public:
	Clock()
	{ 
		dur = 0;
	}
	
	inline void start()
	{
		tmS = time(NULL);
	}
	
	inline void stop()
	{
		dur = time(NULL) - tmS;
	}
	
	inline time_t duration()
	{
		return dur;
	}
	
	inline string stamp()
	{
		// get time
		time_t now = time((time_t*) NULL); 
		struct tm *to = localtime(&now);
		
		// update values
		to->tm_year += 1900;
		to->tm_mon  += 1;
		
		// generate filename
		sprintf(strbuf, "%02d:%02d:%02d", to->tm_hour, to->tm_min, to->tm_sec);
		return string(strbuf) + string(" ");
	}
};


// log file
struct Msg
{
	// write title
	static void title()
	{
#ifdef LG_ON
		cout << PRODUCT_NAME << endl;
		cout << PRODUCT_COPY << endl;
		cout << PRODUCT_AUTHORS << endl;
		next();
#endif
	}
	
	static void line()
	{
#ifdef LG_ON
		cout << endl;
#endif
	}

	static void write(const char* pszFormat, ...)
	{
#ifdef LG_ON
		// buffer and timer
		static char buf[2048];
		static Clock cl;

		// set time stamp
		cout << cl.stamp().c_str();
							
		// write the formated log string to szLine
		va_list argList;
		va_start( argList, pszFormat );
		vsprintf( buf, pszFormat, argList );
		va_end( argList );
		cout << buf << endl;
#endif
	}

	static void rewrite(const char* pszFormat, ...)
	{
#ifdef LG_ON
		// buffer and timer
		static char buf[2048];
		static Clock cl;

		// set time stamp
		cout << cl.stamp().c_str();
							
		// write the formated log string to szLine
		va_list argList;
		va_start( argList, pszFormat );
		vsprintf( buf, pszFormat, argList );
		va_end( argList );
		cout << buf << "                                                                                               \r" << flush;
#endif
	}

	// error handle
	template < typename T >
	static void error (const char* title, T msg)
	{
        	cout << endl << "================================================================================" << endl << endl;
		cout << "Error occured : " << title << endl;
	        cout << endl << "================================================================================" << endl << endl;
		cout << " " << msg;
        	cout << endl << "================================================================================" << endl << endl;
		
		// exit
		exit(0);
	}

	static void next()
	{
#ifdef LG_ON
		cout << endl << "================================================================================" << endl << endl;
#endif
	}
};

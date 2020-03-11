#define		ERROR_WARNING		1
#define		ERROR_CRITICAL		2

#define		IDEFIX_ERROR(ERROR_MESSAGE)		ErrorHandler(ERROR_CRITICAL, ERROR_MESSAGE, std::string(__func__) , __LINE__, std::string(__FILE__))
#define     IDEFIX_WARNING(WARNING_MESSAGE) ErrorHandler(ERROR_WARNING, WARNING_MESSAGE, std::string(__func__) , __LINE__, std::string(__FILE__))

void ErrorHandler (const int,
			  std::stringstream &, 
			  std::string , 
			  const int, 
			  std::string  );

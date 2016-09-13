#include <string.h>
#include <opencv/cv.h>
#include <opencv/cxcore.h>

#define X_VALUE 0
#define Y_VALUE 1

class InputQuadrangle
{
	public:
		double *mApexs;
	public:
		InputQuadrangle(void);
		~InputQuadrangle(void);
		void setFirstPointToRectification(double x,double y);
		void setSecondPointToRectification(double x,double y);
		void setThirdPointToRectification(double x,double y);
		void setFourthPointToRectification(double x,double y);
};


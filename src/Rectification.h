#include "math.h"
#include <string.h>
#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include "InputQuadrangle.h"

#define MAX_WIDTH_OR_HEIGHT 800
#define POINTS 4
#define CORRECTION_VALUE 0.001

class Rectification
{
	public:
		CvMat *mHMatrix;
		double mWidth;
		double mHeight;
		/**********************************new********************/
	private:
		double *newPictureXApex;
		double *newPictureYApex;
		double *mPictureApexXPosition;
		double *mPictureApexYPosition;
		double mOldFirstApexX;
		double mOldFirstApexY;
		CvMat* mFirstApex;
		CvMat* mSecondApex;
		CvMat* mThirdApex;
		CvMat* mFourthApex;
		CvMat* mN2Matrix;
		CvMat* mN3Matrix;
	public:
		Rectification(void);
		~Rectification(void);
		void calculateDimension(InputQuadrangle* quadrangle);
		IplImage* rectify(IplImage* input,IplImage* out);
		/**********************************new********************/
	public:
		IplImage* correctPersperctive(InputQuadrangle* quadrangle,IplImage* input,IplImage* output);
	private:
		void setStartMatrixs(InputQuadrangle* quadrangle);
		CvMat* calculateCenterPoint();
		double calculateFocalLength(double u0, double v0);
		double calculateAspectRatio(double f, double u0, double v0);
		void findProperSizeOfNewImage(double aspectRatio);
		void rewriteApexsToMatrixs();
		IplImage* rectifyImage(IplImage* in, IplImage* out);
		void correctHomographicMatrix(IplImage* inImage,CvMat* invH);
		CvMat* calculateHomographicMatrix();
};


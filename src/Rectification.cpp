#include "Rectification.h"


Rectification::Rectification(void):mHMatrix(NULL),mPictureApexXPosition(NULL),mPictureApexYPosition(NULL),
			newPictureXApex(NULL),newPictureYApex(NULL),mWidth(0),mHeight(0)
{
	mFirstApex = cvCreateMat(3, 1, CV_32FC1);
	mSecondApex = cvCreateMat(3, 1, CV_32FC1);
	mThirdApex = cvCreateMat(3, 1, CV_32FC1);
	mFourthApex = cvCreateMat(3, 1, CV_32FC1);
	mN2Matrix = cvCreateMat(3,1,CV_32FC1);
	mN3Matrix = cvCreateMat(3,1,CV_32FC1);
}


Rectification::~Rectification(void)
{
	if(mHMatrix!=NULL) {
		cvReleaseMat(&mHMatrix);
	}
	if(mPictureApexXPosition!=NULL) {
		delete[] mPictureApexXPosition;
	}
	if(mPictureApexYPosition!=NULL) {
		delete[] mPictureApexYPosition;
	}
	if(newPictureXApex!=NULL) {
		delete[] newPictureXApex;
	}
	if(newPictureYApex!=NULL) {
		delete[] newPictureYApex;
	}
	cvReleaseMat(&mFirstApex);
	cvReleaseMat(&mSecondApex);
	cvReleaseMat(&mThirdApex);
	cvReleaseMat(&mFourthApex);
	cvReleaseMat(&mN2Matrix);
	cvReleaseMat(&mN3Matrix);
}

IplImage* Rectification::correctPersperctive(InputQuadrangle* quadrangle,IplImage* input,IplImage* output) {
	return NULL;
}

void Rectification::setStartMatrixs(InputQuadrangle* quadrangle) {
	/**remember first apex value, it will be needed for correction later*/
	mOldFirstApexX = quadrangle->mApexs[1*2+X_VALUE]-CORRECTION_VALUE;
	mOldFirstApexY = quadrangle->mApexs[1*2+Y_VALUE]-CORRECTION_VALUE;

	/**set apexs values, below apex order*/
	/**1-----------2**********/
	/**-           -**********/
	/**3-----------4**********/
	cvmSet(mSecondApex,0,0,quadrangle->mApexs[2*2+X_VALUE]-mOldFirstApexX);
	cvmSet(mSecondApex,1,0,quadrangle->mApexs[2*2+Y_VALUE]-mOldFirstApexY);
	cvmSet(mSecondApex,2,0,1);
	cvmSet(mThirdApex,0,0,quadrangle->mApexs[0*2+X_VALUE]-mOldFirstApexX);
	cvmSet(mThirdApex,1,0,quadrangle->mApexs[0*2+Y_VALUE]-mOldFirstApexY);
	cvmSet(mThirdApex,2,0,1);
	cvmSet(mFourthApex,0,0,quadrangle->mApexs[3*2+X_VALUE]-CORRECTION_VALUE-mOldFirstApexX);
	cvmSet(mFourthApex,1,0,quadrangle->mApexs[3*2+Y_VALUE]-CORRECTION_VALUE-mOldFirstApexY);
	cvmSet(mFourthApex,2,0,1);
	/**first apex should be at (0,0) point*/
	cvmSet(mFirstApex,0,0,mOldFirstApexX-mOldFirstApexX);
	cvmSet(mFirstApex,1,0,mOldFirstApexY-mOldFirstApexY);
	cvmSet(mFirstApex,2,0,1);
}

CvMat* Rectification::calculateCenterPoint() {
	CvMat *centerMatrix,*addApexs14,*addApexs23,*twoMatrix,*divApexs23,*divApexs14;
	/**allocate memory for matrix*/
	centerMatrix = cvCreateMat(3,1,CV_32FC1);
	addApexs14 = cvCreateMat(3,1,CV_32FC1);
	addApexs23 = cvCreateMat(3,1,CV_32FC1);
	twoMatrix = cvCreateMat(3,1,CV_32FC1);
	divApexs23 = cvCreateMat(3,1,CV_32FC1);
	divApexs14 = cvCreateMat(3,1,CV_32FC1);
	/**fill matrix */
	/* [2]*/
	/* [2]*/
	/* [2]*/
	/*******/
	cvmSet(twoMatrix,0,0,2);
	cvmSet(twoMatrix,1,0,2);
	cvmSet(twoMatrix,2,0,2);

	/**add first apex to fourth*/
	/** [x1] + [x4] = [x1+x4]*/
	/** [y1] + [y4] = [y1+y4]*/
	/** [z1] + [z4] = [z1+z4]*/
	/*****************/
	cvAdd(mFirstApex,mFourthApex,addApexs14);
	/**add second apex to third*/
	/** [x2] + [x3] = [x2+x3]*/
	/** [y2] + [y3] = [y2+y3]*/
	/** [z2] + [z3] = [z2+z3]*/
	/*****************/
	cvAdd(mSecondApex,mThirdApex,addApexs23);
	/**divide previous addition 14 by twoMatrix*/
	cvDiv(addApexs14,twoMatrix,divApexs14);
	/**divide previous addition 23 by twoMatrix*/
	cvDiv(addApexs23,twoMatrix,divApexs23);
	/**find intersection point between two central diagonal points*/
	/**add */
	cvAdd(divApexs14,divApexs23,centerMatrix);
	/**divide by twoMatrix*/
	cvDiv(centerMatrix,twoMatrix,centerMatrix);

	/**release*/
	cvReleaseMat(&divApexs14);
	cvReleaseMat(&divApexs23);
	cvReleaseMat(&twoMatrix);
	cvReleaseMat(&addApexs23);
	cvReleaseMat(&addApexs14);
	return centerMatrix;
}

double Rectification::calculateFocalLength(double u0, double v0) {
	CvMat *cross14,*cross24,*cross34;
	cross14 = cvCreateMat(3, 1, CV_32FC1);
	cross24 = cvCreateMat(3, 1, CV_32FC1);
	cross34 = cvCreateMat(3, 1, CV_32FC1);
	/**calculate cross product of apexs*/
	/** a1 x a4 = cross14*/
	cvCrossProduct(mFirstApex,mFourthApex,cross14);
	/** a2 x a4 = cross24*/
	cvCrossProduct(mSecondApex,mFourthApex,cross24);
	/** a3 x a4 = cross34*/
	cvCrossProduct(mThirdApex,mFourthApex,cross34);

	/**next etap: calculate dot product and divide by another dot product*/
	double k2 = cvDotProduct(cross14,mThirdApex)/cvDotProduct(cross24,mThirdApex);
	double k3 = cvDotProduct(cross14,mSecondApex)/cvDotProduct(cross34,mSecondApex);

	cvReleaseMat(&cross14);
	cvReleaseMat(&cross24);
	cvReleaseMat(&cross34);

	CvMat *mul_k2_apex2,*mul_k3_apex3,*k2_mat,*k3_mat;
	k2_mat = cvCreateMat(3,1,CV_32FC1);
	cvmSet(k2_mat,0,0,k2);
	cvmSet(k2_mat,1,0,k2);
	cvmSet(k2_mat,2,0,k2);
	k3_mat = cvCreateMat(3,1,CV_32FC1);
	cvmSet(k3_mat,0,0,k3);
	cvmSet(k3_mat,1,0,k3);
	cvmSet(k3_mat,2,0,k3);
	mul_k2_apex2 = cvCreateMat(3,1,CV_32FC1);
	mul_k3_apex3 = cvCreateMat(3,1,CV_32FC1);
	cvMul(k2_mat,mSecondApex,mul_k2_apex2);
	cvMul(k3_mat,mThirdApex,mul_k3_apex3);
	cvmSub(mul_k2_apex2,mFirstApex,mN2Matrix);
	cvmSub(mul_k3_apex3,mFirstApex,mN3Matrix);

	cvReleaseMat(&mul_k2_apex2);
	cvReleaseMat(&mul_k3_apex3);
	cvReleaseMat(&k2_mat);
	cvReleaseMat(&k3_mat);

	double f_pow = (-1/(cvmGet(mN2Matrix,2,0)*cvmGet(mN3Matrix,2,0)))*
			((cvmGet(mN2Matrix,0,0)*cvmGet(mN3Matrix,0,0)-(cvmGet(mN2Matrix,2,0)*cvmGet(mN3Matrix,2,0)+cvmGet(mN2Matrix,2,0)*cvmGet(mN3Matrix,0,0))*u0
			+cvmGet(mN2Matrix,2,0)*cvmGet(mN3Matrix,2,0)*u0*u0)+(cvmGet(mN2Matrix,1,0)*cvmGet(mN3Matrix,1,0)-(cvmGet(mN2Matrix,1,0)*cvmGet(mN3Matrix,2,0)+cvmGet(mN2Matrix,2,0)*cvmGet(mN3Matrix,1,0))*v0
			+cvmGet(mN2Matrix,2,0)*cvmGet(mN3Matrix,2,0)*v0*v0));

	double f=750;
	if(f_pow > 0) {
		f = sqrt(f_pow);
	}
	return f;
}

double Rectification::calculateAspectRatio(double f, double u0, double v0) 
{
	/**find aspect ratio (width/height) of new image*/
	CvMat *matrix_A,*n2_trans,*a_trans,*a_inv,*a_trans_inv,*n3_trans,*mul_temp,*mul_temp2,*value,*value2;
	matrix_A = cvCreateMat(3,3,CV_32FC1);
	a_trans = cvCreateMat(3,3,CV_32FC1);
	a_inv = cvCreateMat(3,3,CV_32FC1);
	a_trans_inv = cvCreateMat(3,3,CV_32FC1);
	n2_trans = cvCreateMat(1,3,CV_32FC1);
	n3_trans = cvCreateMat(1,3,CV_32FC1);
	mul_temp = cvCreateMat(1,3,CV_32FC1);
	mul_temp2 = cvCreateMat(1,3,CV_32FC1);
	value = cvCreateMat(1,1,CV_32FC1);
	value2 = cvCreateMat(1,1,CV_32FC1);
	cvmSet(matrix_A,0,0,f);
	cvmSet(matrix_A,0,1,0);
	cvmSet(matrix_A,0,2,u0);
	cvmSet(matrix_A,1,0,0);
	cvmSet(matrix_A,1,1,f);
	cvmSet(matrix_A,1,2,v0);
	cvmSet(matrix_A,2,0,0);
	cvmSet(matrix_A,2,1,0);
	cvmSet(matrix_A,2,2,1);
	cvTranspose(mN2Matrix,n2_trans);
	cvTranspose(mN3Matrix,n3_trans);
	cvTranspose(matrix_A,a_trans);
	cvInvert(matrix_A,a_inv);
	cvInvert(a_trans,a_trans_inv);
	cvmMul(n2_trans,a_trans_inv,mul_temp);
	cvmMul(mul_temp,a_inv,mul_temp2);
	cvmMul(mul_temp2,mN2Matrix,value);
	cvmMul(n3_trans,a_trans_inv,mul_temp);
	cvmMul(mul_temp,a_inv,mul_temp2);
	cvmMul(mul_temp2,mN3Matrix,value2);
	cvReleaseMat(&mul_temp2);
	cvReleaseMat(&mul_temp);
	cvReleaseMat(&n3_trans);
	cvReleaseMat(&a_inv);
	cvReleaseMat(&a_trans);
	cvReleaseMat(&a_trans_inv);
	cvReleaseMat(&n2_trans);
	cvReleaseMat(&matrix_A);
	double aspectRatio_pow = cvmGet(value,0,0)/cvmGet(value2,0,0);
	cvReleaseMat(&value);
	cvReleaseMat(&value2);
	double aspectRatio = sqrt(aspectRatio_pow);
	if(aspectRatio < 0.4) {
		aspectRatio = 0.4;
	}
	if(aspectRatio> 2.5) {
		aspectRatio = 2.5;
	}
	return aspectRatio;
}

void Rectification::findProperSizeOfNewImage(double aspectRatio) 
{
	/**sizes of old quadrangle image*/
	double W1 = sqrt(pow(cvmGet(mSecondApex,0,0)-cvmGet(mFirstApex,0,0),2)+pow(cvmGet(mSecondApex,1,0)-cvmGet(mFirstApex,1,0),2));
	double W2 = sqrt(pow(cvmGet(mFourthApex,0,0)-cvmGet(mThirdApex,0,0),2)+pow(cvmGet(mFourthApex,1,0)-cvmGet(mThirdApex,1,0),2));
	double H1 = sqrt(pow(cvmGet(mThirdApex,0,0)-cvmGet(mFirstApex,0,0),2)+pow(cvmGet(mThirdApex,1,0)-cvmGet(mFirstApex,1,0),2));
	double H2 = sqrt(pow(cvmGet(mFourthApex,0,0)-cvmGet(mSecondApex,0,0),2)+pow(cvmGet(mFourthApex,1,0)-cvmGet(mSecondApex,1,0),2));

	mWidth = W1;
	if(mWidth<W2)
		mWidth = W2;
	mHeight = H1;
	if(mHeight<H2)
		mHeight = H2;

	double oldAspectRatio = mWidth/mHeight;
	if (oldAspectRatio >= aspectRatio) {
		mHeight = mWidth / aspectRatio;
	}
	else {
		mWidth = mHeight * aspectRatio;
	}

	/**not necessary condition*/
	while((mWidth>MAX_WIDTH_OR_HEIGHT)||(mHeight>MAX_WIDTH_OR_HEIGHT)) {
		mWidth = (int)(mWidth/2);
		mHeight = (int)(mHeight/2);
	}
	mWidth = (int)(mWidth);
	mHeight = (int)(mHeight);
}

void Rectification::rewriteApexsToMatrixs() 
{
	mPictureApexXPosition = new double[4];
	mPictureApexXPosition[0] = cvmGet(mFirstApex,0,0);
	mPictureApexXPosition[1] = cvmGet(mSecondApex,0,0);
	mPictureApexXPosition[2] = cvmGet(mThirdApex,0,0);
	mPictureApexXPosition[3] = cvmGet(mFourthApex,0,0);
	mPictureApexYPosition = new double[4];
	mPictureApexYPosition[0] = cvmGet(mFirstApex,1,0);
	mPictureApexYPosition[1] = cvmGet(mSecondApex,1,0);
	mPictureApexYPosition[2] = cvmGet(mThirdApex,1,0);
	mPictureApexYPosition[3] = cvmGet(mFourthApex,1,0);
	newPictureXApex = new double[4];
	newPictureXApex[0] = 0;
	newPictureXApex[1] = mWidth;
	newPictureXApex[2] = 0;
	newPictureXApex[3] = mWidth;
	newPictureYApex = new double[4];
	newPictureYApex[0] = 0;
	newPictureYApex[1] = 0;
	newPictureYApex[2] = mHeight;
	newPictureYApex[3] = mHeight;
}

IplImage* Rectification::rectifyImage(IplImage* inImage, IplImage* outImage) 
{
	double x1 = 0;
	double x2 = 0;
	double x3 = 0;
	double y_position = 0;
	double x_position = 0;
	uchar* data = (uchar *)inImage->imageData;
	uchar* dataOut = (uchar *)outImage->imageData;
	for(int row=0;row<mHeight;row++) {

		for(int col=0;col<mWidth;col++) {

			x1 = cvmGet(mHMatrix,0,0) * col + cvmGet(mHMatrix,0,1) * row + cvmGet(mHMatrix,0,2);
			x2 = cvmGet(mHMatrix,1,0) * col + cvmGet(mHMatrix,1,1) * row + cvmGet(mHMatrix,1,2);
			x3 = cvmGet(mHMatrix,2,0) * col + cvmGet(mHMatrix,2,1) * row + 1;
			y_position = x2/x3 + mOldFirstApexY;

			if(inImage->height < y_position) {
				y_position = (inImage->height-1);
			}

			x_position = x1/x3 + mOldFirstApexX;
			if(inImage->width < x_position) {
				x_position = (inImage->width-1);
			}

			int temp_y = (int)y_position ;
			int temp_x = (int)x_position ;

			if(dataOut!=NULL && data!=NULL) {
				dataOut[row*outImage->widthStep+col*outImage->nChannels] = data[temp_y*inImage->widthStep+temp_x*inImage->nChannels];
				dataOut[row*outImage->widthStep+col*outImage->nChannels+1] = data[temp_y*inImage->widthStep+temp_x*inImage->nChannels+1];
				dataOut[row*outImage->widthStep+col*outImage->nChannels+2] = data[temp_y*inImage->widthStep+temp_x*inImage->nChannels+2];
			}
		}
	}
	return outImage;
}

void Rectification::correctHomographicMatrix(IplImage* inImage,CvMat* invH)
{
	CvMat *hCoeff = cvCreateMat(3,3,CV_32FC1);
	CvMat* multipleMat = cvCreateMat(3,3,CV_32FC1);
	double old_height = inImage->height;
	double y_position = mPictureApexYPosition[3];
	double x_position = mPictureApexXPosition[3];
	for (int i=1;i<10;i++) {
		double x1 = cvmGet(invH,0,0) * (x_position) + cvmGet(invH,0,1) * (y_position) + cvmGet(invH,0,2);
		double x2 = cvmGet(invH,1,0) * (x_position) + cvmGet(invH,1,1) * (y_position) + cvmGet(invH,1,2);
		double x3 = cvmGet(invH,2,0) * (x_position) + cvmGet(invH,2,1) * (y_position) + 1;
		x1 = x1/x3;
		x2 = x2/x3;

		double H_coeff = ((x1/mWidth)+(x2/mHeight))/2 - 0.01;
		for(int coeffRow=0;coeffRow<3;coeffRow++) {
			for(int coeffCol=0;coeffCol<3;coeffCol++) {
				cvmSet(hCoeff,coeffRow,coeffCol,H_coeff);
				cvmSet(multipleMat,coeffRow,coeffCol,cvmGet(mHMatrix,coeffRow,coeffCol));
			}
		}
		cvMul(multipleMat,hCoeff,mHMatrix);
		cvInvert(mHMatrix,invH);
	}

	cvReleaseMat(&multipleMat);
	cvReleaseMat(&hCoeff);
}

CvMat* Rectification::calculateHomographicMatrix() 
{
	CvMat* mmat = cvCreateMat(3,3,CV_32FC1);
	CvMat* a = cvCreateMat(POINTS*2,9,CV_32FC1);
	for(int count=1;count<POINTS+1;count++) {
		cvmSet(a,2*count-2,0,newPictureXApex[count-1]);
		cvmSet(a,2*count-2,1,newPictureYApex[count-1]);
		cvmSet(a,2*count-2,2,1);
		cvmSet(a,2*count-2,3,0);
		cvmSet(a,2*count-2,4,0);
		cvmSet(a,2*count-2,5,0);
		cvmSet(a,2*count-2,6,(-newPictureXApex[count-1]*mPictureApexXPosition[count-1]));
		cvmSet(a,2*count-2,7,(-mPictureApexXPosition[count-1]*newPictureYApex[count-1]));
		cvmSet(a,2*count-2,8,-mPictureApexXPosition[count-1]);
		cvmSet(a,2*count-1,0,0);
		cvmSet(a,2*count-1,1,0);
		cvmSet(a,2*count-1,2,0);
		cvmSet(a,2*count-1,3,newPictureXApex[count-1]);
		cvmSet(a,2*count-1,4,newPictureYApex[count-1]);
		cvmSet(a,2*count-1,5,1);
		cvmSet(a,2*count-1,6,(-newPictureXApex[count-1]*mPictureApexYPosition[count-1]));
		cvmSet(a,2*count-1,7,(-mPictureApexYPosition[count-1]*newPictureYApex[count-1]));
		cvmSet(a,2*count-1,8,-mPictureApexYPosition[count-1]);
	}
	CvMat* U  = cvCreateMat(8,8,CV_32FC1);
	CvMat* D  = cvCreateMat(8,9,CV_32FC1);
	CvMat* V  = cvCreateMat(9,9,CV_32FC1);
	CvMat* V22 = cvCreateMat(3,3,CV_32FC1);
	mHMatrix = cvCreateMat(3,3,CV_32FC1);
	CvMat* invH = cvCreateMat(3,3,CV_32FC1);
	cvSVD(a, D, U, V, CV_SVD_U_T|CV_SVD_V_T);

	for(int a=0;a<3;a++) {
		for(int b=0;b<3;b++) {
			cvmSet(mmat,a,b,cvmGet(V,8,a*3+b));
			cvmSet(V22,a,b,(1/cvmGet(V,8,4)));
		}
	}

	cvMul(mmat,V22,mHMatrix);
	cvInvert(mHMatrix,invH);
	cvReleaseMat(&U);
	cvReleaseMat(&D);
	cvReleaseMat(&V);
	cvReleaseMat(&V22);
	cvReleaseMat(&a);
	cvReleaseMat(&mmat);
	return invH;
}

void Rectification::calculateDimension(InputQuadrangle* quadrangle) {
	setStartMatrixs(quadrangle);
	CvMat* imageCenter = calculateCenterPoint();
	double u0 = cvmGet(imageCenter,1,0);
	double v0 = cvmGet(imageCenter,2,0);
	cvReleaseMat(&imageCenter);
	double f = calculateFocalLength(u0,v0);
	double aspectRatio = calculateAspectRatio(f,u0,v0);
	findProperSizeOfNewImage(aspectRatio);
	rewriteApexsToMatrixs();
}

IplImage* Rectification::rectify(IplImage* inImage,IplImage* outImage) {
	CvMat* invH = calculateHomographicMatrix();
	correctHomographicMatrix(inImage,invH);
	rectifyImage(inImage,outImage);
	cvReleaseMat(&invH);
	return outImage;
}

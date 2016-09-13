#include "InputQuadrangle.h"


InputQuadrangle::InputQuadrangle(void): mApexs(NULL)
{
	mApexs = new double[4*2];
	/*********picture apexs***************/
	/***-------------------------------***/
	/***-    1                       2-***/
	/***-                             -***/
	/***-    4                       3-***/
	/***-------------------------------***/
	/*************************************/
	//tablica
	mApexs[0] = 65;
	mApexs[1] = 62;
	mApexs[2] = 397;
	mApexs[3] = 36;
	mApexs[4] = 378;
	mApexs[5] = 334;
	mApexs[6] = 61;
	mApexs[7] = 286;
}


InputQuadrangle::~InputQuadrangle(void)
{
	if(mApexs!=NULL)
		delete[] mApexs;
}

void InputQuadrangle::setFirstPointToRectification(double x,double y) {
	mApexs[0] = x;
	mApexs[1] = y;
}

void InputQuadrangle::setSecondPointToRectification(double x,double y) {
	mApexs[2] = x;
	mApexs[3] = y;
}
	
void InputQuadrangle::setThirdPointToRectification(double x,double y) {
	mApexs[4] = x;
	mApexs[5] = y;
}
	
void InputQuadrangle::setFourthPointToRectification(double x,double y) {
	mApexs[6] = x;
	mApexs[7] = y;
}
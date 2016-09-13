#include <iostream>
//#include <conio.h>
#include <string.h>
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"
#include "ml.h"
#include "Rectification.h"

InputQuadrangle* quadrangle;
int globalMouseCounter = 0;
void mouseHandler(int event, int x, int y, int flags, void *param);

int main() {
	std::string pathSource = "note.jpg";
	IplImage* image = cvLoadImage(pathSource.c_str());
	if(image!=NULL) {
		quadrangle = new InputQuadrangle();
		cvNamedWindow("image",0);
		cvShowImage("image",image);
		cvSetMouseCallback( "image", mouseHandler, NULL );
		while(globalMouseCounter<4) {
			cvWaitKey();
		}
		Rectification *rectification = new Rectification();                             
		rectification->calculateDimension(quadrangle);
		IplImage* output = cvCreateImage(cvSize(rectification->mWidth, rectification->mHeight),IPL_DEPTH_8U,3);
		output = rectification->rectify(image,output);
		cvNamedWindow("image2");
		cvShowImage("image2",output);
		cvWaitKey();
		cvReleaseImage(&output);
		cvReleaseImage(&image);
		cvDestroyWindow("image");
		cvDestroyWindow("image2");
	}
}

void mouseHandler(int event, int x, int y, int flags, void *param)

{
    switch(event) {
        case CV_EVENT_LBUTTONDOWN: {
            fprintf(stdout, "Left button down (%d, %d).\n", x, y);
			if(globalMouseCounter==0) {
				quadrangle->setFirstPointToRectification(x,y);
			}
			else if(globalMouseCounter==1) {
				quadrangle->setSecondPointToRectification(x,y);
			}
			else if(globalMouseCounter==2) {
				quadrangle->setThirdPointToRectification(x,y);
			}
			else if(globalMouseCounter==3) {
				quadrangle->setFourthPointToRectification(x,y);
			}
			globalMouseCounter++;
            break;
		}
    }
}

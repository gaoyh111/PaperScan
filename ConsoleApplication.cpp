
#include "stdafx.h"
#include "highgui.h"
#include <opencv2/opencv.hpp>
#include<time.h>
#include<stdio.h>
using namespace std;
using namespace cv;
void fun_binary(IplImage *gray, float n, IplImage *pDst);
void Contours(IplImage *binary,int** a);
void fun_PerspectiveTransform(int **b, IplImage* gmi, IplImage* pDist);
void fun_colorShadow(IplImage* src, IplImage* pdst);
void colorRestoration(IplImage* merged, IplImage* dst, IplImage* pDst);
void histogramEqualization(IplImage* src, IplImage* pDst);

int main()  {
	clock_t start, med, finish, tt1, tt2, tt3, tt4, tt5;
	double t1, t2, t3, t4, t5, t6, t7, t8, totaltime;//计时用的													
	start = clock();

	//char path[100] = "C:/Users/hasee/Desktop/C++Projects/test2.jpg";
	char path[100] = "C:/Users/hasee/Desktop/C++Projects/source4-48.jpg";
	IplImage *resize = NULL;
	IplImage *img = NULL, *pic = NULL; CvSize sz;
	pic = cvLoadImage(path, 1);//读取图像  //0.78s


	if (pic->width > 2000)//判断尺寸大小,大的话默认按照3042*4032,尺寸缩小3倍到1008,1344
	{
		resize = cvCreateImage(cvSize(pic->width / 3, pic->height / 3), pic->depth, pic->nChannels);
		cvResize(pic, resize, CV_INTER_CUBIC);
	}
	else
	{
		resize = cvCreateImage(cvSize(pic->width, pic->height), pic->depth, pic->nChannels);
		cvResize(pic, resize, CV_INTER_CUBIC);//resize是尺寸适中的图片
	}

	if (resize->width > resize->height)//读取图片的bug,某些牌子的手机拍的照片读入长宽会颠倒,判断并矫正
	{
		img = cvCreateImage(cvSize(resize->height, resize->width), resize->depth, resize->nChannels);
		cvResize(resize, img, CV_INTER_CUBIC);//长宽互换
		cvTranspose(resize, img);                                                        //0.15s
		cvFlip(img, NULL, 1);//进行旋转
		resize = cvCreateImage(cvSize(resize->height, resize->width), resize->depth, resize->nChannels);
		cvCopy(img, resize);//得到正常的图片
	}

	sz.width = resize->width*0.5;//尺寸再缩小2倍,二值化操作减少计算量
	sz.height = resize->height*0.5;
	IplImage* resizeSmaller = cvCreateImage(sz, resize->depth, resize->nChannels);//resizeSmaller尺寸很小,二值图计算量减少很多
	cvResize(resize, resizeSmaller, CV_INTER_CUBIC);
	IplImage* gray = cvCreateImage(cvGetSize(resizeSmaller), resizeSmaller->depth, 1);
	cvCvtColor(resizeSmaller, gray, CV_BGR2GRAY);//resizeSmaller的灰度图
	tt1 = clock();//计时.在二值图开始之前,是图相读取,准备的耗时.

	fun_binary(gray, 0.93, gray);
	/*
	fun_binary(IplImage *gray, float n, IplImage *pDst);
	函数实现的功能：对图相进行二值化处理
	输入：
	IplImage* gray：灰度图
	float n:阈值(大小滤波的比值)
	输出：
	IplImage* dst：二值图
	*/


	//对二值图寻找轮廓
	int **b = NULL;//创建一个二维数组b,分配[5][2]个空间,用来存放4个角点,长和宽.
	b = (int **)malloc(5 * sizeof(int *));
	for (int i = 0; i<5; i++)
		b[i] = (int *)malloc(2 * sizeof(int));
	Contours(gray, b);
	/*
	void Contours(IplImage *binary, int** a);
	函数实现的功能：对图相进行轮廓查找,筛选出目标轮廓,获得其角点,其最小外接矩形的长和宽
	输入：
	IplImage* binary：二值图
	输出：
	int** a：二维数组
	*/

	if (b[3][0] == 0 && b[3][1] == 0)//没有角点,则为初始值0,0 
	{
		cout << "未检测到角点!";
	}
	tt2 = clock();

	int w, h;
	w = 2 * b[4][1]; //最小外接矩形的宽;
	h = 2 * b[4][0];//高
	IplImage* pTransPic = 0;
	pTransPic = cvCreateImage(cvSize(w, h), IPL_DEPTH_8U, 3);//创建最小矩形的宽高的图像,接收透视变换处理好的图

	fun_PerspectiveTransform(b, resize, pTransPic);
	/*
	void fun_PerspectiveTransform(int **b, IplImage* gmi, IplImage* pDist);
	函数实现的功能：对图像进行透视变换
	输入：
	int **b:包含角点,最小外接矩形长宽信息的二维数组
	IplImage* gmi:尺寸合适的原图
	输出：
	IplImage* pDist：透视变换后的图
	*/
	tt3 = clock();//透视变换的耗时
	sz.width = w;//尺寸指定为最小外接矩形的长宽,也可以指定为定值###846*1146###
	sz.height = h;



	//↓修改的地方.主函数末尾加上cvReleaseImage(&pDest);释放pDest
	//之前变量名ROI改成了pDest,然后创建一个感兴趣ROI,是pDest边距减少7像素的图ROI,之后均对ROI进行处理.
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	IplImage* pDest = 0;//pDest接收四角定位后的图片                                                                   //
	pDest = cvCreateImage(sz, resize->depth, resize->nChannels);                                                      //
	cvResize(pTransPic, pDest, CV_INTER_CUBIC);//修改尺寸                                                             //
	med = clock();// 四角定位的总耗时                                                                                 //
	                                                                                                                  //
	                    //具体裁剪多少像素,可自行修改                                                                 //
	CvSize size = cvSize(pDest->width , pDest->height);//从坐标(7,7)到(width - 15, pDest->height - 15)ROI区域//
	cvSetImageROI(pDest, cvRect(0,0, size.width, size.height));//设置源图像ROI                                       //
	IplImage* ROI = cvCreateImage(size, pDest->depth, pDest->nChannels);//创建目标图像(之后处理的图像                 //
	cvResize(pDest, ROI);                                                                                             //
	//cvCopy(ROI, pDest); //复制图像                                                                                  //
	//cvResetImageROI(pDest);//源图像用完后，清空ROI(用不上)                                                          //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	IplImage *temp = cvCreateImage(cvGetSize(ROI), IPL_DEPTH_8U, 3);
	cvCopy(ROI, temp);
	IplImage* bgRemove = cvCreateImage(cvGetSize(ROI), IPL_DEPTH_8U, 3);//bgRemove是背景去除后的图片
	fun_colorShadow(temp, bgRemove);//去除背景
	/*
	void fun_colorShadow(IplImage* src, IplImage* pdst);
	函数实现的功能：对图像进行背景去除
	输入：
	IplImage* src:原图
	输出：
	IplImage* pDst:去除背景后的图
	*/
	tt4 = clock();


	IplImage* colorRecover = cvCreateImage(cvGetSize(ROI), IPL_DEPTH_8U, 3);
	cvCopy(bgRemove, colorRecover);
	colorRestoration(bgRemove, ROI, colorRecover);//还原色彩
	/*
	void colorRestoration(IplImage* bgRemove, IplImage* ROI, IplImage* pDst);
	函数实现的功能：对图像进行色彩还原
	输入：
	IplImage* bgRemove:去除背景的图片
	IplImage* ROI:包含背景的图片
	输出：
	IplImage* pDst:还原色彩的图片
	*/
	tt5 = clock();


	IplImage* enhanced = cvCreateImage(cvGetSize(ROI), IPL_DEPTH_8U, 3);;
	histogramEqualization(colorRecover, enhanced);//直方图均衡//enhanced是颜色增强后的图片
	/*
	void histogramEqualization(IplImage* src, IplImage* pDst);
	函数实现的功能：对图像进行色彩还原
	输入：
	IplImage* src:原图
	输出：
	IplImage* pDst:彩色增强后的图 方法:直方图均衡
	*/
	finish = clock();
	t3 = (double)(tt1 - start) / CLOCKS_PER_SEC;//计算时间
	t4 = (double)(tt2 - tt1) / CLOCKS_PER_SEC;
	t5 = (double)(tt3 - tt2) / CLOCKS_PER_SEC;
	t6 = (double)(tt4 - tt3) / CLOCKS_PER_SEC;
	t7 = (double)(tt5 - tt4) / CLOCKS_PER_SEC;
	t8 = (double)(finish - tt5) / CLOCKS_PER_SEC;
	t1 = (double)(med - start) / CLOCKS_PER_SEC;
	t2 = (double)(finish - med) / CLOCKS_PER_SEC;
	totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "\n读取图片" << t3 << "秒！";
	cout << "\n二值化" << t4 << "秒！";
	cout << "\n透视变换" << t5 << "秒！";
	cout << "\n角点检测运行时间为" << t1 << "秒！";
	cout << "\n背景去除" << t6 << "秒！";
	cout << "\n彩色还原" << t7 << "秒！";
	cout << "\n彩色增强" << t8 << "秒！";
	cout << "\n背景增强运行时间为" << t2 << "秒！";
	cout << "\n此程序的运行时间 为" << totaltime << "秒！" << endl;
	
	cvShowImage("gray", enhanced);
	cvWaitKey(0);

	cvReleaseImage(&resize);
	cvReleaseImage(&resizeSmaller);
	cvReleaseImage(&pic);
	cvReleaseImage(&img);
	cvReleaseImage(&gray);
	cvReleaseImage(&ROI);
	cvReleaseImage(&pTransPic);
	cvReleaseImage(&temp);
	cvReleaseImage(&pDest);
	cvReleaseImage(&bgRemove);
	cvReleaseImage(&colorRecover);
	cvReleaseImage(&enhanced);
	delete[] b;
		
	
	return 0;
}
void fun_binary(IplImage *gray, float n, IplImage *pDst)//二值化
{
	IplImage* bg = cvCreateImage(cvGetSize(gray), IPL_DEPTH_8U, 1);
	IplImage* top = cvCreateImage(cvGetSize(gray), IPL_DEPTH_8U, 1);
	cvSmooth(gray, top, CV_GAUSSIAN, 3, 3, 0, 0);//进行小滤波处理,每个像素值接近前景
	cvSmooth(gray, bg, CV_GAUSSIAN, 21, 21, 0, 0);//进行大滤波处理,平均了像素,每个像素值更接近背景
	cvScale(bg, bg, 0.93, 0.0);//bg=bg*0.93    
	cvCmp(top, bg, pDst, CV_CMP_GT);//top与bg*0.93比较.也就是top/bg和阈值0.93的比较
	cvReleaseImage(&bg);
	cvReleaseImage(&top);
}
void Contours(IplImage *pBinary,int** a)//查找轮廓,返回角点
{

	a[0][0] = 0; a[0][1] = 0;//初始化二维数组
	a[1][0] = 0; a[1][1] = 0;
	a[2][0] = 0; a[2][1] = 0;
	a[3][0] = 0; a[3][1] = 0;
	a[4][0] = 0; a[4][1] = 0;

	CvSeq *pContour = NULL;
	CvSeq *pConInner = NULL;
	CvMemStorage *pStorage = NULL;
	pStorage = cvCreateMemStorage(0);
	CvMemStorage* point = cvCreateMemStorage();
	CvSeq* mcont = NULL; double area = 10000000; int seqSize = 0;
	CvBox2D box; int w, h; int temp = 0;double dConArea; int count = 0;
	IplImage * copy = cvCreateImage(cvGetSize(pBinary), IPL_DEPTH_8U, 1);
	cvCopy(pBinary, copy);//finContours函数会改动传入的参数,用临时变量代替
	cvFindContours(copy, pStorage, &pContour, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);//查找轮廓,pContour是指向第一个轮廓的指针,借此用循环遍历所有轮廓
	for (CvSeq * contour = pContour; contour != nullptr; contour = contour->h_next)
	{
		dConArea = fabs(cvContourArea(contour,CV_WHOLE_SEQ));//计算轮廓面积
		box=cvMinAreaRect2(contour);//轮廓的最小外接矩形
		w = box.size.height;
		h = box.size.width;
		if (dConArea > 0.196 *pBinary->width*pBinary->height && dConArea < 0.738 *pBinary->width*pBinary->height && w*h>0.5*dConArea&&w*h<1.5*dConArea)//根据面积大小初步筛选
		{
			double e = 0.1*cvArcLength(contour);
			mcont=cvApproxPoly(contour, sizeof(CvContour), point, CV_POLY_APPROX_DP, e, 0);//金丝多边形函数,把轮廓近似成多边形
			//cvDrawContours(resultImage,mcont, cvScalar(255, 0, 0), cvScalar(0, 0, 255), 0, 2, 8);
			seqSize = mcont->total;//角点个数
			if (seqSize == 4 && area > dConArea &&w*h>0.5*dConArea && w*h < 1.5*dConArea)//判断近似多边形的角点个数,筛选角点
			{
				area = dConArea;
				CvSeqBlock* firstElement = mcont->first;
				CvPoint* curp = (CvPoint*)firstElement->data;
				for (int j = 0; j < seqSize; j++)
				{
					a[j][0] = curp->x;//把角点坐标存放到数组a中
					a[j][1] = curp->y;
					//cout << a[j][0]<<","<<a[j][1] << endl;
					++curp;
				}
				if (w<h)
				{
					temp = w;
					w = h;
					h = temp;
				}
				a[4][0] = w; a[4][1] = h;//把长,宽存放到数组a中
				
			}
			
		}
	}
	cvReleaseMemStorage(&pStorage);
	cvReleaseMemStorage(&point);
	cvReleaseImage(&copy);
}
void fun_PerspectiveTransform(int **b, IplImage* gmi, IplImage* pDst)//透视变换, 先确定角点的具体位置,再进行透视变换
{
	
	cvZero(pDst);
	int c[4][2];
	int w, h;
	w = 2 * b[4][1]; //最小矩形的宽;
	h = 2 * b[4][0];//高
	int min = 20000000;
	int max = 1;
	int l, a1, a2, a3, a4;
	for (int i = 0; i < 4; i++)//根据角点(x,y)中x+y的值可以确定出距离远点最近和最远的两个点
	{
		l = b[i][0] + b[i][1];
		if (min>l)
		{
			min = l;
			a1 = i;
		}
		if (max < l)
		{
			max = l;
			a4 = i;
		}
	}
	c[0][0] = 2*b[a1][0]; c[0][1] = 2*b[a1][1];
	c[3][0] = 2*b[a4][0]; c[3][1] = 2*b[a4][1];
	int hmax = 1;
	int hmin = 20000000;
	for (int i = 0; i < 4; i++)//确定另外两个点的位置
	{
		if (i != a1&&i != a4)
		{
			if (hmax < b[i][1])
			{
				hmax = b[i][1];
				a3 = i;

			}
			if (hmin > b[i][1])
			{
				hmin = b[i][1];
				a2 = i;

			}
		}
	}
	c[1][0] = 2*b[a2][0]; c[1][1] = 2*b[a2][1];
	c[2][0] = 2*b[a3][0]; c[2][1] = 2*b[a3][1];
	CvMat *warp_matrix = cvCreateMat(3, 3, CV_32FC1);
	//透视变换
	CvPoint2D32f srcQuad[4], dstQuad[4];//透视变换的参数,srcQuad存放角点,dstQuad存放变换后的相应位置.
	srcQuad[0].x = c[0][0];
	srcQuad[0].y = c[0][1];
	srcQuad[1].x = c[1][0];
	srcQuad[1].y = c[1][1];
	srcQuad[2].x = c[2][0];
	srcQuad[2].y = c[2][1];
	srcQuad[3].x = c[3][0];
	srcQuad[3].y = c[3][1];
	dstQuad[0].x = 0;
	dstQuad[0].y = 0;
	dstQuad[1].x = w;
	dstQuad[1].y = 0;
	dstQuad[2].x = 0;
	dstQuad[2].y = h;
	dstQuad[3].x = w;
	dstQuad[3].y = h;
	
	
	cvGetPerspectiveTransform(srcQuad, dstQuad, warp_matrix);
	
	cvWarpPerspective(gmi, pDst, warp_matrix);
	cvReleaseMat(&warp_matrix);
}
void fun_colorShadow(IplImage* src, IplImage* pdst)//粗略背景去除
{
	IplImage *B, *G, *R;//RGB三个通道
	IplImage *tempB, *tempG, *tempR;//RGB三个通道的复制
	CvSize Size1 = cvGetSize(src);
	B = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	tempB = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	G = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	tempG = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	R = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	tempR = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	cvSplit(src, B, G, R, 0);
	cvCopy(B, tempB);
	cvCopy(G, tempG);
	cvCopy(R, tempR);
	IplImage* maxB = cvCreateImage(Size1, IPL_DEPTH_8U, 1);//对每个通道进行大小滤波,提取主要轮廓
	IplImage* minB = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	cvSmooth(B, minB, CV_GAUSSIAN, 3, 3, 0, 0);
	cvSmooth(B, maxB, CV_GAUSSIAN, 31, 31, 0, 0);
	IplImage* maxG = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	IplImage* minG = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	cvSmooth(G, minG, CV_GAUSSIAN, 3, 3, 0, 0);
	cvSmooth(G, maxG, CV_GAUSSIAN, 31, 31, 0, 0);
	IplImage* maxR = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	IplImage* minR = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	cvSmooth(R, minR, CV_GAUSSIAN, 3, 3, 0, 0);
	cvSmooth(R, maxR, CV_GAUSSIAN, 31, 31, 0, 0);
	IplImage* A1 = cvCreateImage(Size1, IPL_DEPTH_8U, 1);//做abs(G-B)/G,通过阈值区分背景和前景
	cvAbsDiff(G, B, A1);
	IplImage* A2 = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	cvAbsDiff(G, R, A2);
	for (int i = 0; i < src->height; i++)
	{

		char *pA1 = A1->imageData + i * A1->widthStep;
		char *pA2 = A2->imageData + i * A2->widthStep;
		char *pB = B->imageData + i * B->widthStep;
		char *pG = G->imageData + i * G->widthStep;
		char *pR = R->imageData + i * R->widthStep;
		for (int j = 0; j < src->width; j++)//概要:大小滤波得到主要轮廓,彩色和背景rgb值的差异,区分出背景.
		{//通过if语句区分背景和前景:大小滤波比值大于0.96,abs(tempG-tempB)/tempG<0.16,abs(tempG-tempR)/tempG<0.12,tempG>110
			if (((uchar *)(minB->imageData + i*minB->widthStep))[j]>0.96*((uchar *)(maxB->imageData + i*maxB->widthStep))[j] && pA1[j]<0.16*((uchar *)(tempG->imageData + i*tempG->widthStep))[j] && pA2[j]<0.12*((uchar *)(tempG->imageData + i*tempG->widthStep))[j] && ((uchar *)(tempG->imageData + i*tempG->widthStep))[j]>110)
			{
				pB[j] = 255;
			}
			if (((uchar *)(minG->imageData + i*minG->widthStep))[j]>0.96*((uchar *)(maxG->imageData + i*maxG->widthStep))[j] && pA1[j]<0.16*((uchar *)(tempG->imageData + i*tempG->widthStep))[j] && pA2[j]<0.12*((uchar *)(tempG->imageData + i*tempG->widthStep))[j] && ((uchar *)(tempG->imageData + i*tempG->widthStep))[j]>110)
			{
				pG[j] = 255;
			}
			if (((uchar *)(minR->imageData + i*minR->widthStep))[j]>0.96*((uchar *)(maxR->imageData + i*maxR->widthStep))[j] && pA1[j]<0.16*((uchar *)(tempG->imageData + i*tempG->widthStep))[j] && pA2[j]<0.12*((uchar *)(tempG->imageData + i*tempG->widthStep))[j] && ((uchar *)(tempG->imageData + i*tempG->widthStep))[j]>110)
			{
				pR[j] = 255;
			}
				
		}
	}



	cvMerge(B, G, R, 0, pdst);
	cvReleaseImage(&B);
	cvReleaseImage(&G);
	cvReleaseImage(&R);
	cvReleaseImage(&tempB);
	cvReleaseImage(&tempG);
	cvReleaseImage(&tempR);
	cvReleaseImage(&maxB);
	cvReleaseImage(&minB);
	cvReleaseImage(&maxG);
	cvReleaseImage(&minG);
	cvReleaseImage(&maxR);
	cvReleaseImage(&minR);
	cvReleaseImage(&A1);
	cvReleaseImage(&A2);

	//cvShowImage("gray", dst);
	//cvWaitKey(0);
}
void colorRestoration(IplImage* bgRemove, IplImage* ROI, IplImage* pDst)//还原色彩.查找大轮廓,将其填充到纯黑的二值图中,得到大轮廓的白板,根据这个还原大轮廓的色彩.
{
	CvSize Size1 = cvGetSize(ROI);
	IplImage* gray = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	IplImage* binary = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	IplImage* temp = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	IplImage* finalBin = cvCreateImage(Size1, IPL_DEPTH_8U, 1);//纯黑的二值图
	cvZero(finalBin);
	IplImage* Bin = cvCreateImage(Size1, IPL_DEPTH_8U, 1);//纯黑的二值图
	cvCvtColor(bgRemove, gray, CV_BGR2GRAY);
	cvThreshold(gray, binary, 254, 255, CV_THRESH_BINARY);
	cvCopy(binary, temp);
	
	CvSeq *pContour = NULL;
	CvSeq *pConInner = NULL;
	CvMemStorage *pStorage = NULL;
	pStorage = cvCreateMemStorage(0);
	CvMemStorage* point = cvCreateMemStorage();
	CvSeq* mcont = NULL;
	int w = Size1.width, h = Size1.height;
	int areaS = 0; double dConArea = 0;
	cvFindContours(temp, pStorage, &pContour, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE);//查找轮廓
	for (CvSeq * contour = pContour; contour != nullptr; contour = contour->h_next)
	{
		cvZero(Bin);
		dConArea = fabs(cvContourArea(contour, CV_WHOLE_SEQ));
		if (dConArea > 0.1 * w * h && dConArea < 0.7 * w * h)
		{
			cvDrawContours(Bin, contour, cvScalar(255), cvScalar(255), 0, -1, 8);
			for (int i = 0; i < ROI->height; i++)
			{
				char *pB = Bin->imageData + i * Bin->widthStep;
				char *pb = binary->imageData + i * binary->widthStep;
				for (int j = 0; j < ROI->width; j++)
				{
					if (((uchar *)(Bin->imageData + i*Bin->widthStep))[j] == 255 && ((uchar *)(binary->imageData + i*binary->widthStep))[j] == 255)//计算去除背景后的大轮廓内的白色像素的数量
					{
						areaS++;
					}
				}
			}
			cout << areaS;
			if (areaS < 0.35*dConArea)//白色像素少于一定比例,认为应该还原.
			{
				cvDrawContours(finalBin, contour, cvScalar(255), cvScalar(255), 0, -1, 8);//在finalBin中涂白.
			}	
		}
		
	}



	for (int i = 0; i < ROI->height; i++)
	{
		char *ptrDst = ROI->imageData + i * ROI->widthStep;
		char *pMerged = pDst->imageData + i * pDst->widthStep;
		char *pFb = finalBin->imageData + i * finalBin->widthStep;
		for (int j = 0; j < ROI->width; j++)
		{
			if (pFb[j] != 0)//根据finalBin进行彩色还原
			{
				pMerged[3 * j] = ptrDst[3 * j];
				pMerged[3 * j + 1] = ptrDst[3 * j + 1];
				pMerged[3 * j + 2] = ptrDst[3 * j + 2];
			}
		}
	}


	cvReleaseImage(&gray);
	cvReleaseImage(&binary);
	cvReleaseImage(&temp);
	cvReleaseImage(&finalBin);
	cvReleaseImage(&Bin);
	cvReleaseMemStorage(&pStorage);
	cvReleaseMemStorage(&point);


}
void histogramEqualization(IplImage* src, IplImage* pDst)//直方图均衡,对rgb三个通道做直方图处理,将三个直方图整合到一个,并计算累积分布和出现概率.获得查找表.根据查找表更新rgb值
{
	IplImage *B, *G, *R;
	CvSize Size1 = cvGetSize(src);
	B = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	G = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	R = cvCreateImage(Size1, IPL_DEPTH_8U, 1);
	cvSplit(src, B, G, R, 0);
	int hist_size = 256; float range[] = {5,250 }; float* ranges[] = { range };
	CvHistogram* r_hist = cvCreateHist(1, &hist_size, CV_HIST_ARRAY, ranges, 1);//直方图的初始化
	CvHistogram* g_hist = cvCreateHist(1, &hist_size, CV_HIST_ARRAY, ranges, 1);
	CvHistogram* b_hist = cvCreateHist(1, &hist_size, CV_HIST_ARRAY, ranges, 1);


	cvCalcHist(&B, b_hist, 0, 0);//计算直方图
	cvCalcHist(&G, g_hist, 0, 0);
	cvCalcHist(&R, r_hist, 0, 0);
#define cvQueryHistValue_1D(hist,idx0)((float)cvGetReal1D((hist)->bins, (idx0)))//获取直方图的数值的函数
	int r[256] = { 0 }, g[256] = { 0 }, b[256] = { 0 };
	int s[256] = { 0 };
	int width = src->width;
	int height = src->height;
	int sum = width * height;       
	int i, j;

	//分别统计直方图的RGB分布

	b[0] = cvQueryHistValue_1D(b_hist, 0);
	g[0] = cvQueryHistValue_1D(g_hist, 0);
	r[0] = cvQueryHistValue_1D(r_hist, 0);
	s[0] = b[0] + r[0] + g[0];
	for (i = 1; i < 256; i++)
	{
		b[i] = cvQueryHistValue_1D(b_hist, i) + b[i - 1];
		g[i] = cvQueryHistValue_1D(g_hist, i) + g[i - 1];
		r[i] = cvQueryHistValue_1D(r_hist, i) + r[i - 1];
		s[i] = b[i] + g[i] + r[i];//计算三个直方图的累积分布并且整合到s数组,s数组相当于是彩图的直方图的累积分布.
	}
	

	for (int i = 1; i <256; i++)
	{
		s[i] = s[i] * 255 / (s[255]+1);//均衡化.
	}
	//cvLUT(src, s, dst);
	//归一化直方图
	for (i = 0; i<height; i++)
	for (j = 0; j<width; j++)
	{
		((uchar*)(pDst->imageData + i*pDst->widthStep))[j*pDst->nChannels + 0] = s[((uchar*)(src->imageData + i*src->widthStep))[j*src->nChannels + 0]];//根据s数组更新rgb的数值
		((uchar*)(pDst->imageData + i*pDst->widthStep))[j*pDst->nChannels + 1] = s[((uchar*)(src->imageData + i*src->widthStep))[j*src->nChannels + 1]];
		((uchar*)(pDst->imageData + i*pDst->widthStep))[j*pDst->nChannels + 2] = s[((uchar*)(src->imageData + i*src->widthStep))[j*src->nChannels + 2]];
	}
	//https://blog.csdn.net/tengfei461807914/article/details/76952276
	//直方图均衡的算法参考的是上面的博客
	cvReleaseImage(&B);
	cvReleaseImage(&G);
	cvReleaseImage(&R);
	cvReleaseHist(&r_hist);
	cvReleaseHist(&g_hist);
	cvReleaseHist(&b_hist);
}

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

typedef  unsigned char imgpix;
#define	 PIXVAL	  255
#define  DIMENSION	8
#define  SEEDLEN	4		

static _inline double abs_doub(double x)
{
	return (x<0) ? (-x) : x;
}

static _inline int round(double x)
{
	if (x<0)
		return (int)(x-0.5);
	else
		return (int)(x+0.5);
}

/*
宏定义
*/
#define SIFT_MIN(x, y) ((x<y) ? x : y)

#define interp_hist_peak(l, c, r) (0.5*((l)-(r))/((l)+(r)-2.0*(c)))
#define dsqrt(x, y) sqrt((x*x)+(y*y))
#define datan(y, x) atan2(y, x)
#define fsqrt(sigma, sigman) sqrt(((sigma)*(sigma))-((sigman)*(sigman)))
#define fpow(a, b) pow(a, b)
#define fdsqrt(x) sqrt(x);

//#define floor(x) ((int)(x>0 ? x : (x-1)))

#define MAX_IMG_VAL					255
#define SIFT_GAUSSPYR_INTVL			3
#define SIFT_DOGPYR_LEVEL			(SIFT_GAUSSPYR_LEVEL-1)
#define SIFT_MAX_GAUSS_SMOOTH_SIG	6.4			//6.4*wid==gauss窗口半径
#define SIFT_MAX_INTERP_STEPS		5			//精确定位最大循环次数
#define SIFT_EDGE_RATIO				10.0
#define SIFT_ORI_HIST_BIN_NUM		36
#define SIFT_DES_HIST_BIN_NUM		8
#define SIFT_HIST_WIDTH				4
#define SIFT_DES_HIST_FCT			3			//决定描述子一个直方图所占图像宽度
#define SIFT_ORI_SIG				1.5			//计算梯度方向时，平滑窗口所用尺度的倍数
#define SIFT_ORI_RADIUS				(3.0*SIFT_ORI_SIG)	//region半径的倍数
#define SIFT_PI						3.1415926
#define SIFT_PI2					6.28318531
#define SIFT_SMOOTH_HIST_LEVEL		2			//平滑方向直方图的次数
#define SIFT_HIST_PEAK_RATIO		0.8			//直方图主方向设定比率阈值
#define SIFT_DESC_LEN				128			//4*4*8

/*
結構體
*/
typedef struct
{
	int		width;
	int		height;
	imgpix	**imgdata;
}ImgInfo;


#endif
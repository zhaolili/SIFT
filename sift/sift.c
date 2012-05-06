/************************************************************************
* Author:	Zhao Lili <Beihang University, "zhao86.scholar@gmail.com">
* Version:	V3.3
* TIme:		Aug.,2010	  
************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "matrix.h"
#include "mem.h"
#include "global.h"

#define SIFT_API _declspec(dllexport)
#include "sift.h"

/*
sift

infile:		保留参数
imgdata:	图像数据
width:		图像宽度
height:		图像高度

//information below can be input through UI
wid:		gauss卷积模板半径为gwid*sigma
db_img:		是否double图像
sigma0:		基本层图像的平滑因子		
sigman:		假设原始图像已存在的平滑因子
contr_thr:	低对比度阈值=contr_thr/scales
sift_border_dist: 检测边界距离
ill_norm_thresh:  歸一化閾值
top_oct_reso:     最后组图像分辨率

float wid=2, int db_img=0, float sigma0 =1.6f,\
float sigman=0.5f, float contr_thr=0.09f, int sift_border_dist=5, float ill_norm_thresh=0.2f
*/
FeatureVect* _sift_main(const char* infile, imgpix	*imgData, int width, int height, float wid, int db_img, float sigma0,\
						float sigman, float contr_thr, int sift_border_dist, float ill_norm_thresh, int top_oct_reso)
{
	int octave_num;								//组数
	int w, h;
	int scales = SIFT_GAUSSPYR_INTVL;
	FeatureVect *featVect;

	ImgInfo *imgInf;
	float	****gaussimg;		//存放gauss滤波图像,[octaves][scales][height][width]
	float	****dogimg;			//dog图像
	float	****magimg;
	float	****ortimg;
	
	//获取图像数据
	if (db_img)
	{
		imginfo_fill_dbl(imgData, &imgInf, width, height);
		w = width<<1;
		h = height<<1;
		sigman *= 2;
	}		
	else
	{
		imginfo_fill(imgData, &imgInf, width, height);
		w = width;
		h = height;
	}		
	
	//计算尺度空间组数,最顶层图像的分辨率最小为top_oct_reso
	octave_num = (int)(log(SIFT_MIN(w, h))/log(2)) - (int)(log(top_oct_reso)/log(2))+1;

	//为gauss滤波图像和DOG图像
	gauss_mem_alloc(octave_num, scales+3, imgInf, &gaussimg, &dogimg);
	
	//梯度幅度和方向分配空间
	mag_ort_mem_alloc(&magimg, &ortimg, octave_num, scales, w, h);
	
	//计算得到gauss金字塔空间
	gauss_space(wid, octave_num, scales, sigma0, sigman, gaussimg, imgInf);
	
	//计算得到dog空间
	dog_space(octave_num, scales, dogimg, gaussimg, w, h);

	//得到梯度幅度和方向
	cal_grad_mag_ori(magimg, ortimg, octave_num, scales, w, h, gaussimg);

	//检测dog空间极值点
	featVect = scale_space_extrema(sift_border_dist, dogimg, octave_num, scales, w, h, contr_thr);

	//计算极值点的平滑因子
	cal_feat_sigma(featVect, sigma0, scales);

	//如果第一层的base图像被double了
	if (db_img)
		adjust_doubled_img(featVect);

	//计算关键点主方向
	cal_feat_oris(featVect, magimg, ortimg, imgInf);

	//计算描述子
	cal_descriptor(featVect, magimg, ortimg, imgInf, SIFT_HIST_WIDTH, SIFT_DES_HIST_BIN_NUM, ill_norm_thresh);

	//gauss濾波圖像和dog圖像釋放空間
	gauss_mem_free(octave_num, scales, gaussimg, dogimg);

	//释放梯度幅度和方向空间
	mag_ort_mem_free(magimg, ortimg, octave_num, scales);

	//釋放原始圖像（或double后的）的圖像空間
	imginfo_release(imgInf);

	////釋放特徵點鏈錶
	//free_feature_nodes(featVect);
	return featVect;
}

/*
釋放特徵點鏈錶

featVect: 待釋放鏈錶
*/
void free_feature_nodes(FeatureVect *featVect)
{
	FeatureNode *prev;
	FeatureNode *tail = featVect->cur_node;

	while(tail!=NULL)
	{
		free(tail->feat);
		prev = tail->prev;
		tail->prev = NULL;
		tail->next = NULL;
		free(tail);
		tail = prev;
	}
	featVect->cur_node = NULL;
	featVect->first = NULL;
	free(featVect);
}

/*
计算高斯图像梯度幅度和方向

mag:  梯度空间
ort:  梯度方向空间
octaves: 组数
scales:  实际有用的层数,3
w:	  base组的图像宽度
h:    base组的图像高度
gaussimg: gauss图像
*/
int cal_grad_mag_ori(float ****magimg, float ****ortimg, int octaves, int scales, int w, int h, float ****gaussimg)
{
	int o, s, i, j;
	float dx, dy;
	w = w<<1;
	h = h<<1;
	for (o=0; o<octaves; o++)
	{
		w = w>>1;
		h = h>>1;
		for (s=1; s<=scales; s++)
		{
			for (j=1; j<h-1; j++)	//边除外
			{
				for (i=1; i<w-1; i++)
				{
					dx = gaussimg[o][s][j][i+1]-gaussimg[o][s][j][i-1];
					dy = gaussimg[o][s][j-1][i]-gaussimg[o][s][j+1][i];
					magimg[o][s-1][j][i] = sqrt(dx*dx+dy*dy);
					if ((dx==0.0)&&(dy==0.0))
						ortimg[o][s-1][j][i] = 0.0f;
					else
					{
						ortimg[o][s-1][j][i] = atan2(dy, dx);
						//ortimg[o][s-1][j][i] = (ortimg[o][s-1][j][i]<0)?(ortimg[o][s-1][j][i]+SIFT_PI2):ortimg[o][s-1][j][i];
					}						
				}
			}
		}
		
	}
	return 0;
}

/*
为梯度和方向分配空间

mag:  梯度空间
ort:  梯度方向空间
octaves: 组数
scales:  实际有用的层数,3
w:	  base组的图像宽度
h:    base组的图像高度
*/
int mag_ort_mem_alloc(float *****mag, float *****ort, int octaves, int scales, int w, int h)
{
	int j;
	int mem = 0;

	w = w<<1;
	h = h<<1;

	*mag = (float ****)calloc(octaves, sizeof(float ***));
	*ort = (float ****)calloc(octaves, sizeof(float ***));

	for (j=0; j<octaves; j++)
	{
		w = w>>1;
		h = h>>1;
		mem += mem3D_alloc_float(&(*mag)[j], scales, h, w);
		mem += mem3D_alloc_float(&(*ort)[j], scales, h, w);
	}

	return mem;
}

/*
释放梯度和方向空间

mag:  要释放的梯度空间
ort:  要释放的梯度方向空间
octaves: 组数
scales:  实际有用的层数,3
*/
void mag_ort_mem_free(float ****mag, float ****ort, int octaves, int scales)
{
	mem4D_free_float(mag, octaves, scales);
	mem4D_free_float(ort, octaves, scales);
}


/*
计算描述子.计算结构保存在Feature结构体中的descriptor数组中

featVect:	存储极值点链表信息
magimg:    梯度值图像
ortimg：   梯度方向
imginf:    原始图像信息（或加倍后的）
d:		   直方图数组宽度，默认为4
binnum:	   每个直方图的方向数，默认为8
*/
void cal_descriptor(FeatureVect *featVect, float ****magimg, float ****ortimg, ImgInfo *imginf, int d, int binnum, float thresh)
{
	int i;
	float ***hist;		//直方图，[r][c][octave]
	FeatureNode *node = featVect->first;
	Feature *feat;

	mem3D_alloc_float(&hist, d, d, binnum);
	while(node != NULL)
	{
		feat = node->feat;
		
		//计算得到hist
		cal_desc_hist(hist, magimg, ortimg, imginf, node, d, binnum);
		
		//3D->1D
		desc_1Dto3D(node, hist, d, binnum);
		
		//將descriptor歸一化
		desc_norm(feat->descriptor, SIFT_DESC_LEN);
		
		//去除高梯度值5
		for (i=0; i<SIFT_DESC_LEN; i++)
			if (feat->descriptor[i]>thresh)
				feat->descriptor[i] = thresh;
		
		//再次歸一化
		desc_norm(feat->descriptor, SIFT_DESC_LEN);
		node = node->next;
	}

	mem3D_free_float(hist, d);
}
/*
歸一化描述子,將descriptor歸一化到單位長度

descriptor: 待歸一化的描述子
len:        描述子長度 默認為128
*/
void desc_norm(float *descriptor, int len)
{
	int i;
	float sum=0.0f;
	for (i=0; i<len; i++)
		sum += descriptor[i]*descriptor[i];

	sum = fdsqrt(sum);
	for (i=0; i<len; i++)
		descriptor[i] /= sum;
}
/*
將3維直方圖裝換位一維，存儲在feature中

node:	當前關鍵點結點
hist:	要轉換的直方圖
d:		直方圖寬度
binnum:	直方圖中bin的個數
*/
void desc_1Dto3D(FeatureNode *node, float ***hist, int d, int binnum)
{
	int j, i,t;
	Feature *feat = node->feat;

	for (j=0; j<d; j++)
	{
		for (i=0; i<d; i++)
		{
			t = (j*d+i)*binnum;
			memcpy(feat->descriptor+t, hist[j][i], binnum*sizeof(float));
		}
	}
}

/*
计算得到descriptor的直方图

hist:		存储结果[hang][lie][bin]		
magimg:     梯度值图像
ortimg:     方向
imginf:     原始图像信息（或加倍后的）
Node:		当前结点
d:		    直方图数组宽度，默认为4
binnum:	    每个直方图的方向数，默认为8
*/
void cal_desc_hist(float ***hist, float ****magimg, float ****ortimg, ImgInfo *imginf, FeatureNode *Node, int d, int binnum)
{
	int			i,j, times, radius;
	int			octave, s, x, y, w, h;
	int			r0, c0, o0, r, c, o, rb, cb, ob;
	double		cos_t, sin_t, r_rot, c_rot, rbin, cbin, obin;
	double		hist_width, mag, sigma, keyori,rot_ori, exp_denom, win, bins_per_rad;
	double		d_r, d_c, d_o, v_r, v_c, v_o;
	Feature *feat = Node->feat;

	octave = feat->octave;
	s = feat->scale;
	x = feat->x;
	y = feat->y;
	times = (int)fpow(2.0, octave);
	h = imginf->height/times;
	w = imginf->width/times;
	sigma		= feat->sig_octave;			//精确到亚尺度
	keyori		= feat->orient;
	hist_width	= SIFT_DES_HIST_FCT*sigma;
	radius		= (int)(sqrt(2.0)*hist_width*(d+1.0)*0.5+0.5);	//描述子所占图像半径
	exp_denom	= d*d*0.5;					//2*0.5*d*0.5*d
	bins_per_rad= binnum/SIFT_PI2;

#define fast_floor(x) ((int)( x - ((x>=0)?0:1) ))

	//reset hist
	for (j=0; j<d; j++)
		for (i=0; i<d; i++)
			memset(hist[j][i], 0, sizeof(float)*binnum);

	//计算主方向角度的三角函数值
	cos_t = cos(keyori);
	sin_t = sin(keyori);

	for (j=-radius; j<=radius; j++)
	{
		if ((y-j>h-2)||(y-j<1))
			continue;
		for (i=-radius; i<=radius; i++)
		{
			if ((x+i>w-2)||(x+i<1))   //边上
				continue;

			//计算旋转后的坐标，以关键点为圆心,以直方图为单位
			c_rot = (j*cos_t-i*sin_t)/hist_width;
			r_rot = (j*sin_t+i*cos_t)/hist_width;
			rbin  = r_rot + d/2 -0.5;
			cbin  = c_rot + d/2 -0.5;

			//梯度和角度
			mag = magimg[octave][s-1][y-j][x+i];
			rot_ori = ortimg[octave][s-1][y-j][x+i];

			//在范围内
			if (rbin>-1.0 && rbin<d && cbin>-1.0 && cbin<d)
			{
				//计算当前像素旋转后角度
				rot_ori = rot_ori-keyori;
				while(rot_ori<0.0)
					rot_ori = rot_ori+SIFT_PI2;
				while(rot_ori>=SIFT_PI2)
					rot_ori = rot_ori-SIFT_PI2;

				//计算当前像素处于直方图中的位置
				obin = rot_ori*bins_per_rad;

				//高斯权
				win = exp(-(c_rot*c_rot+r_rot*r_rot)/exp_denom);

				//將當前像素的梯度值vote到直方圖的8個bin中
				r0 = floor(rbin);
				c0 = floor(cbin);
				o0 = floor(obin);
				d_r = rbin-r0;
				d_c = cbin-c0;
				d_o = obin-o0;

				//trilinear interpolation
				for (r=0; r<=1; r++)
				{
					rb = r0 + r;
					if (rb>=0 && rb<d)
					{
						v_r =  win*mag*((r==0) ? 1.0-d_r : d_r);
						for (c=0; c<=1; c++)
						{
							cb = c0+c;
							if (cb>=0 && cb<d)
							{
								v_c = v_r*((c==0) ? 1.0-d_c : d_c);
								for (o=0; o<=1; o++)
								{
									ob = (o0+o)%binnum;
									if (ob>=binnum)
									{
										printf("");
									}
									v_o = v_c*((o==0) ? 1.0-d_o : d_o);
									hist[rb][cb][ob] += v_o;
								}
							}
						}						
					}
				}//for
			}
		}
	}
}


/*
计算关键点主方向

featVect:	存储极值点链表信息
magimg:     梯度幅度值
ortimg:     梯度方向值
imginf:    原始图像信息（或加倍后的）
*/
void cal_feat_oris(FeatureVect *featVect, float ****magimg, float ****ortimg, ImgInfo *imginf)
{
	int i;
	int	x, y, w, h, times;
	int	octave, scale;
	int 		rad;		//计算关键点主方向时的region半径
	float		sig;		//在当前层的精确尺度
	Feature		*feat;
	FeatureNode *node;
	double *hist = (double *)calloc(SIFT_ORI_HIST_BIN_NUM, sizeof(double));
	double omax;		//直方图最大幅值

	w = imginf->width;
	h = imginf->height;
	node = featVect->first;
	while(node!=NULL)
	{
		feat	= node->feat;
		sig		= SIFT_ORI_SIG*feat->sig_octave;
		rad		= (int)(SIFT_ORI_RADIUS*feat->sig_octave + 0.5);

		//计算当前关键点的梯度直方图
		x = feat->x;
		y = feat->y;
		octave = feat->octave;
		scale  = feat->scale;
		times  = (int)fpow(2.0, octave);
		w = w/times;
		h = h/times;
		cal_hist(magimg[octave][scale-1], ortimg[octave][scale-1], y, x, rad, sig, hist, w, h);

		//平滑方向直方图
		for (i=0; i<SIFT_SMOOTH_HIST_LEVEL; i++)
			smooth_hist(hist, SIFT_ORI_HIST_BIN_NUM);

		//计算直方图中最大幅度值
		omax = domin_mag(hist, SIFT_ORI_HIST_BIN_NUM);

		//计算关键点主方向，更新链表 返回下一个结点的指针
		node = add_dominated_ori(hist, omax*SIFT_HIST_PEAK_RATIO, SIFT_ORI_HIST_BIN_NUM, node, featVect);
	}
	free(hist);
}

/*
计算关键点主方向，并将幅值>=omax的方向也作为关键点主方向（将其作为一个新关键点插入列表中）

hist:		直方图数组
mag_thr:	幅度阈值，omax*SIFT_HIST_PEAK_RATIO
bins:		直方图中的bin数
featNode:	当前关键点
featVect:	存储极值点链表信息
返回值：	下一个结点
*/
FeatureNode* add_dominated_ori(double *hist, double mag_thr, int bins, FeatureNode *featNode, FeatureVect *featVect)
{
	int i, signal = 0;
	int l, r;
	double nb, ori;

	Feature		*feat;
	FeatureNode *newNode, *nextNode, *retNode = featNode->next;

	for (i=0; i<bins; i++)
	{
		l = (i==0)?(bins-1):i-1;
		r = (i+1==bins) ? 0 : i+1;
		if((hist[i]>hist[l])&&(hist[i]>hist[r])&&(hist[i]>=mag_thr))
		{
			//计算插值后的方向
			nb = i + interp_hist_peak(hist[l], hist[i], hist[r]);
			nb = (nb<0) ? (bins+nb) : (nb>=bins ? (nb-bins) : nb);
			ori = ((SIFT_PI2*nb)/bins) - SIFT_PI;
			signal++;
			if (signal==1)
			{
				feat = featNode->feat;
				feat->orient = ori;
			}
			else
			{
				feat = new_feature();
				feat_copy(feat, featNode->feat);
				feat->orient = ori;
				newNode = (FeatureNode *)calloc(1, sizeof(FeatureNode));
				newNode->feat = feat;
				//将新结点插入链表中
				nextNode = featNode->next;
				if (nextNode!=NULL)
					nextNode->prev = newNode;
				newNode->next = nextNode;
				newNode->prev = featNode;
				featNode->next = newNode;	
				featVect->to_extr_num++;		//关键点个数增加
			}			
		}
	}
	return retNode;
}

/*
拷贝特征点

dest: 指向目标feature结构
src:  指向源feature结构
*/
void feat_copy(Feature *dest, Feature *src)
{
	dest->octave = src->octave;
	dest->ox = src->ox;
	dest->oy = src->oy;
	dest->scale = src->scale;
	dest->sig_octave = src->sig_octave;
	dest->sig_space = src->sig_space;
	dest->sub_sca = src->sub_sca;
	dest->x = src->x;
	dest->y = src->y;
}

/*
计算dominated 关键点幅度

hist: 直方图
bins: 直方图中的bin数
返回值： 最大幅度值
*/
double domin_mag(double *hist, int bins)
{   
	int i;
	int binmax;
	double omax;

	binmax = 0; 
	omax = hist[0];
	for (i=1; i<bins; i++)
	{
		if (hist[i]>omax)
		{
			omax	= hist[i];
			binmax	= i;
		}
	}
	return omax;
}

/*
平滑方向直方图

hist: 待平滑的直方图
bins: 直方图中的bin数
*/
void smooth_hist(double *hist, int bins)
{
	int i;
	double prev, next, tmp;

	prev = hist[bins-1];
	for (i=0; i<bins; i++)
	{
		tmp = hist[i];
		next = (i+1 == bins) ? hist[0] : hist[i+1];
		hist[i] = 0.25*prev + 0.5*tmp + 0.25*next;
		prev = tmp;
	}
}

/*
计算当前关键点的梯度直方图

magimg: 梯度值
ortimg: 方向
y:     关键点纵坐标
x:     关键点横坐标
rad:   region半径
sigma: region平滑尺度
hist:  直方图
w:     当前关键点所在octave的图像宽度
h:     图像高度
*/
void cal_hist(float **magimg, float **ortimg, int y, int x, int rad, float sigma, double *hist, int w, int h)
{
	int i, j;
	int bin;
	double mag, ori, exp_denom, wei;

	memset(hist, 0, sizeof(float)*SIFT_ORI_HIST_BIN_NUM);

	exp_denom = 2.0*sigma*sigma;
	for (j=-rad; j<=rad; j++)
	{
		if ((y-j<1)||(y-j>h-2))
			continue;
		for (i=-rad; i<=rad; i++)
		{
			if ((x+i<1)||(x+i>w-2))		//不在范围内
				continue;
			mag = magimg[y-j][x+i];
			ori = ortimg[y-j][x+i];
			wei = exp(-(j*j+i*i)/exp_denom);
			bin = round(SIFT_ORI_HIST_BIN_NUM*(ori+SIFT_PI)/SIFT_PI2);
			bin = (bin<SIFT_ORI_HIST_BIN_NUM) ? bin : 0;		//[0,10).....[340, 350)[350, 360)
			hist[bin] += wei*mag;
		}
	}
}


/*
若第一层的base图像被double了，修改feature信息中的全局坐标和平滑因子

featVect:	存储极值点链表信息
*/
void adjust_doubled_img(FeatureVect *featVect)
{
	Feature		*feat;
	FeatureNode *node;

	node = featVect->first;
	while(node!=NULL)
	{
		feat = node->feat;
		feat->ox /= 2.0;
		feat->oy /= 2.0;
		feat->sig_space /= 2.0;

		node = node->next;
	}
}

/*
计算所有极值点的尺度

featVect:	存储极值点链表信息
sigma:      base层的尺度
scales：    层数
*/
void cal_feat_sigma(FeatureVect *featVect, float sigma, int scales)
{
	int octave;
	float scale; //精确到亚位置
	FeatureNode *node = featVect->first;
	Feature *feat;

	while(node!=NULL)
	{
		feat	= node->feat;
		octave	= feat->octave;
		scale	= feat->scale+feat->sub_sca;

		feat->sig_octave = sigma*fpow(2.0, scale/scales);
		feat->sig_space = feat->sig_octave*fpow(2.0, octave);

		node = node->next;
	}
}

/*
得到所有极值点（精确定位、去除低对比度及边缘点后的）

bd_dist:  边界点，在此区域内不进行极值检测
dogimg:   dog空间
octaves:  尺度空间组数
scales:   层数
basWid:   base层图像宽度
basHei:   base层图像高度
返回值:   关键点链表信息
*/
FeatureVect* scale_space_extrema(int bd_dist, float ****dogimg, int octaves, int scales, int basWid, int basHei, float contr_thr)
{
	int i, j;
	int o, s;
	int w = basWid<<1;
	int h = basHei<<1;
	float ctr_thresh = contr_thr/scales;		//低对比度阈值
	float ini_thresh = ctr_thresh*0.5f;			//预筛选值
	Feature		*nfeat;
	FeatureNode *featnode;
	FeatureVect *featVect = (FeatureVect *)calloc(1, sizeof(FeatureVect));
	featVect->to_extr_num = 0;
	featVect->first = NULL;
	featVect->cur_node = NULL;

	for (o=0; o<octaves; o++)
	{
		w = w>>1; h = h>>1;
		for (s=1; s<=scales; s++)
		{
			for (j=bd_dist; j<h-bd_dist; j++)
			{
				for (i=bd_dist; i<w-bd_dist; i++)
				{
					//预筛选
					if (dogimg[o][s][j][i]<ini_thresh)
						continue;
					//如果当前位置为极值点
					if (is_local_extrema(dogimg[o], s, j, i))
					{
						//精确定位，并去除低对比度及边界点
						nfeat = acc_localization(w, h, scales, bd_dist, dogimg[o], o, s, j, i, ctr_thresh, SIFT_EDGE_RATIO, SIFT_MAX_INTERP_STEPS);
						if (nfeat)
						{
							featnode = (FeatureNode *)calloc(1, sizeof(FeatureNode));
							featVect->to_extr_num++;			//极值点数++
							featnode->feat = nfeat;				//当前结点中的极值点信息
							featnode->next = NULL;				//指向链表中的下一个结点
							featnode->prev = featVect->cur_node;//指向链表中的前一个结点
							if (featVect->to_extr_num==1)		//若当前结点是第一个结点
								featVect->first = featnode;
							else
								featVect->cur_node->next = featnode;
							featVect->cur_node = featnode;		//将这个极值点作为链表中的当前结点
						}
					}
				}
			}
		}
	}

	return featVect;
}

/*
获得精确定位、去除对低对比度和边界点后的点

w:          图像宽度
h:          图像高度
scales:		层数
bd_dist:    边界点，在此区域内不进行极值检测
dog:		dog图像[scale][y][x]
octave:     所在组
scale:		当前点所在层
y:			待比较图像纵坐标
x:			待比较图像横坐标
ctr_thresh:	对比度门限
ratio:     比率，去除边界点时用，默认为10.0
max_steps:	精确定位最大步长
*/
Feature* acc_localization(int w, int h, int scales, int bd_dist, float ***dog, int octave, int scale, int y, int x, float ctr_thresh, float ratio, int max_steps)
{
	int		cnt = 0, times;
	int		ny, nx, ns;		//新位置
	double	yi, xi, si;		//偏移位置
	double	dy, dx, ds;		//一阶偏导	
	double  dExtr;			//极值
	Feature *feat;

	ns = scale;
	ny = y;
	nx = x;
	while(cnt<max_steps)
	{
		//计算偏移值
		cal_interp_offset(dog, ns, ny, nx, &yi, &xi, &si, &dy, &dx, &ds);
		//找到精确点，跳出循环
		if ((abs_doub(yi)<0.5)&&(abs_doub(xi)<0.5)&&(abs_doub(si)<0.5))
			break;
		nx += round(xi);
		ny += round(yi);
		ns += round(si);

		//新位置是否超出边界范围
		if ((nx<bd_dist)||(ny<bd_dist)||(nx>=w-bd_dist)||(ny>=h-bd_dist)||(ns<1)||(ns>scales))
			return NULL;

		cnt++;
	}
	if (cnt==max_steps)
		return NULL;

	//去除低对比度
	dExtr = dx*xi + dy*yi + ds*si;
	dExtr = dog[ns][ny][nx]+0.5*(dExtr);
	if (dExtr<ctr_thresh)
		return NULL;


	//去除边界点
	if(is_on_edge(dog[scale], ny, nx, ratio))
		return NULL;

	feat = new_feature();
	feat->octave = octave;
	feat->scale = ns;
	feat->y = ny;
	feat->x = nx;
	feat->sub_sca = si;
	times = (int)fpow(2.0, octave);
	feat->ox = (nx+xi)*times;
	feat->oy = (ny+yi)*times;

	return feat;
}

/*
新建一个feature空间
*/
Feature* new_feature()
{
	Feature *feat = (Feature *)calloc(1, sizeof(Feature));
	return feat;
}


/*
是否为边界点

dog:   dog空间[y][x]
y:     纵坐标
x:     横坐标
ratio：最大特征与最小特征值的比率
*/
int is_on_edge(float **dog, int y, int x, float ratio)
{
	double dxx, dxy, dyy;
	double tr, det;
	float  val= dog[y][x];

	dxx = dog[y][x+1]-val - (val-dog[y][x-1]);
	dyy = dog[y+1][x]-val - (val-dog[y-1][x]);
	dxy = ((dog[y+1][x+1]-dog[y+1][x-1])-(dog[y-1][x+1]-dog[y-1][x-1]))/4.0;

	tr	= dxx + dyy;
	det	= dxx*dyy-dxy*dxy;
	if (det<=0)
		return 1;  //在边上
	ratio = ((ratio+1)*(ratio+1))/ratio;
	if(((tr*tr)/det)<ratio)
		return 0;  //低于比率，不在边上的可能性大

	return 1;
}

/*
计算偏移值

dog:		dog图像[scale][y][x]
scale:		当前点所在层
y:			待比较图像纵坐标
x:			待比较图像横坐标
yi:         y偏移
xi:			x偏移
si:			尺度偏移
dy，dx，ds: 偏导
*/
void cal_interp_offset(float ***dog, int scale, int y, int x, double *yi, double *xi, double *si, double *pdy, double *pdx, double *pds)
{
	//使用double类型提高精度
	double dx, dy, ds;		//一阶偏导
	double dxx, dyy, dss, dxy, dxs, dys;//二阶偏导
	double **m;				//矩阵元素
	float  val = dog[scale][y][x];

	dx = (dog[scale][y][x+1]-dog[scale][y][x-1])/2.0;
	dy = (dog[scale][y+1][x]-dog[scale][y-1][x])/2.0;
	ds = (dog[scale+1][y][x]-dog[scale-1][y][x])/2.0;

	dxx = dog[scale][y][x+1]-val - (val-dog[scale][y][x-1]);
	dyy = dog[scale][y+1][x]-val - (val-dog[scale][y-1][x]);
	dss = dog[scale+1][y][x]-val - (val-dog[scale-1][y][x]);
	dxy = ((dog[scale][y+1][x+1]-dog[scale][y+1][x-1])-(dog[scale][y-1][x+1]-dog[scale][y-1][x-1]))/4.0;
	dxs = ((dog[scale+1][y][x+1]-dog[scale+1][y][x-1])-(dog[scale-1][y][x+1]-dog[scale-1][y][x-1]))/4.0;
	dys = ((dog[scale+1][y+1][x]-dog[scale+1][y-1][x])-(dog[scale-1][y-1][x]-dog[scale-1][y-1][x]))/4.0;

	//建立矩阵
	mem2D_alloc_double(&m, 3, 3);
	m[0][0] = dxx;
	m[1][1] = dyy;
	m[2][2] = dss;
	m[0][1] = m[1][0] = dxy;
	m[0][2] = m[2][0] = dxs;
	m[1][2] = m[2][1] = dys;
	if (inverse(m, 3)==1)
	{
		//矩阵乘法
		*xi = m[0][0]*dx + m[0][1]*dy + m[0][2]*ds;
		*yi = m[1][0]*dx + m[1][1]*dy + m[1][2]*ds;
		*si = m[2][0]*dx + m[2][1]*dy + m[2][2]*ds;
	}
	else
		*xi = *yi = *si = 0.0;

	*pdy = dy;
	*pdx = dx;
	*pds = ds;

	mem2D_free_double(m);
}

/*
判断当前点是否为极值点

dog:  dog图像[scale][y][x]
scale:当前点所在层
y:	  待比较图像纵坐标
x:    待比较图像横坐标
*/
int is_local_extrema(float ***dog, int scale, int y, int x)
{
	int s, signal=0;
	int i, j;
	float val = dog[scale][y][x];

	if (val>0)
	{
		for (s=-1; s<=1; s++)
			for (j=-1; j<=1; j++)
				for (i=-1; i<=1; i++)
					if (dog[s+scale][y+j][x+i]>val)
						return 0;
	}
	else
	{
		for (s=-1; s<=1; s++)
			for (j=-1; j<=1; j++)
				for (i=-1; i<=1; i++)
					if (dog[s+scale][y+j][x+i]<val)
						return 0;
	}
	return 1;
}

/*
计算得到dog空间

octaves:  尺度空间组数
scales:	  gauss空间每组层数, +3之前
dogimg:   dog图像
gaussimg: gauss图像
baseWid:  base层的宽度
baseHei:  base层的高度
*/
void dog_space(int octaves, int scales, float ****dogimg, float ****gaussimg, int baseWid, int baseHei)
{
	int o, s;
	int i, j;
	int w = baseWid<<1;
	int h = baseHei<<1;

	for (o=0; o<octaves; o++)
	{
		w = w>>1;
		h = h>>1;
		for (s=0; s<scales+2; s++)		//dog空间比gauss空间少一层
		{
			for (j=0; j<h; j++)
			{
				for (i=0; i<w; i++)
				{
					dogimg[o][s][j][i] = gaussimg[o][s+1][j][i]-gaussimg[o][s][j][i];
				}
			}
		}
	}
}

/*
计算得到gauss尺度空间

wid:      计算gauss平滑模板的窗口大小，wid*sigma,wid默认值为2
octaves:  尺度空间组数
scales:	  gauss空间每组层数, +3之前
sigma0:   最底层的gauss平滑尺度，默认为1.6
sigman:   假设的原始图像已被平滑的尺度, 程序中设为0.5或1.0（doubled）
gaussimg: 高斯平滑图像
imginfo:  原始图像数据（或加倍后的）
*/
void gauss_space(float wid, int octaves, int scales, float sigma0, float sigman, float ****gaussimg, ImgInfo *imginf)
{
	int i,j,o,s,t;
	int gwid;		//gauss窗口半径大小
	int w	= imginf->width;
	int h	= imginf->height;
	int nw, nh;			//边界填充后的宽高
	imgpix  **orgdata = imginf->imgdata;
	float *sigma = (float *)calloc(scales+3, sizeof(float));
	float sigma_prev;		//上一层的尺度
	float ki;
	int max_gwid = (int)(wid*SIFT_MAX_GAUSS_SMOOTH_SIG); 
	float **gsKernel;		//存放高斯核
	float **smoothdata;		//存放滤波数据
	max_gwid = max_gwid<<1;
	max_gwid++;  //最大窗口大小，不是半径
	mem2D_alloc_float(&gsKernel, max_gwid, max_gwid);
	mem2D_alloc_float(&smoothdata, max_gwid+h, max_gwid+w);

	//归一化图像像素值
	for (j=0; j<h; j++)
		for (i=0; i<w; i++)
			gaussimg[0][0][j][i] = ((float)orgdata[j][i])/MAX_IMG_VAL;
			
	//计算各层的平滑尺度
	sigma[0] = fsqrt(sigma0, sigman);
	ki = fpow(2.0, 1.0/scales);
	for (s=1; s<(scales+3); s++)
	{
		sigma_prev = fpow(ki, s-1)*sigma0;
		sigma[s] = sigma_prev*ki;
		sigma[s] = fsqrt(sigma[s], sigma_prev);
	}

	//计算gauss尺度空间
	w = w<<1;
	h = h<<1;
	for (o=0; o<octaves; o++)
	{
		w = w>>1; h = h>>1;
		for (s=0; s<scales+3; s++)
		{
			gwid = (int)(sigma[s]*wid+0.5f);
			gauss_kernel(sigma[s], gwid, gsKernel);
			nw = w + gwid*2;
			nh = h + gwid*2;
			if ((s>=scales)||((s==0)&&(o==0))||(o==0)) 
			{
				if (s!=0)
					t = s-1;
				else
					t = s;
				//将待滤波数据拷贝到临时数据区，进行边界填充
				for (j=0; j<h; j++)
				{
					memcpy(smoothdata[j+gwid]+gwid, gaussimg[o][t][j], w*sizeof(float));
					for (i=0; i<gwid; i++)
					{
						smoothdata[j+gwid][i]	= gaussimg[o][t][j][0];
						smoothdata[j+gwid][nw-i-1]	= gaussimg[o][t][j][w-1];
					}
				}
				for (j=0; j<gwid; j++)
				{
					memcpy(smoothdata[j], smoothdata[gwid], nw*sizeof(float));
					memcpy(smoothdata[nh-j-1], smoothdata[nh-gwid-1], nw*sizeof(float));
				}

				//滤波
				for (j=0; j<h; j++)
				{
					for (i=0; i<w; i++)
					{
						gaussimg[o][s][j][i] = gauss_filter(smoothdata, gsKernel, i+gwid, j+gwid, gwid);
					}
				}
			}
			else  //其他组的第一到scale层
			{
				down_sample(gaussimg[o][s], gaussimg[o-1][scales+s], w, h);	//其他组的第一层的数据直接由上一组的尺度为2sigma的平滑图像下采样得到
			}
		}
	}
	
	mem2D_free_float(gsKernel);
	mem2D_free_float(smoothdata);
}

/*
下采样（最邻近插值法）

dest:		采样后的数据
src:		源数据
destWidth:	目的图像宽
destHeight：目的图像高
*/
void down_sample(float **dest, float **src, int destWidth, int destHeight)
{
	int i, j;
	int i1, j1;

	for(j=0; j<destHeight; j++)
	{
		for (i=0; i<destWidth; i++)
		{
			i1 = i<<1;
			j1 = j<<1;
			dest[j][i] = src[j1][i1];
		}
	}
}

/*
计算高斯核

sig:  平滑尺度
gwid: 平滑半径
gw：  高斯模板
*/
void gauss_kernel(float sig, int gwid, float **gw)
{
	int i, j;
	float sig2 = 2*sig*sig;

	for (j=-gwid; j<=gwid; j++)
	{
		for (i=-gwid; i<=gwid; i++)
			gw[j+gwid][i+gwid] = exp(-(i*i+j*j)/sig2);
	}
}


/*
gauss filter

data: 待滤波数据
gw：  高斯模板
c:    列坐标，i
r:    行坐标, j
w:    滤波半径
*/
float gauss_filter(float **data, float **gw, int c, int r, int w)
{
	int		i, j;
	float	res = 0.0, sum=0.0;
	float	tmp; 

	c = c-w;
	r = r-w;
	w = w+w;
	for (j=0; j<=w; j++)
	{
		for (i=0; i<=w; i++)
		{
			tmp = gw[j][i];
			res += tmp * (data[r+j][c+i]);
			sum += tmp;
		}
	}
	return res/sum;
}


/*
为尺度空间分配内存

octave: 尺度空间组数
scales: 每组层数
imginf: 原始图像（或加倍后的原始图像）
gaussimg: gauss平滑图像[octaves][scales][height][width]
dogimg:	  dog图像
*/
int gauss_mem_alloc(int octaves, int scales, ImgInfo *imginf, float *****gaussimg, float *****dogimg)
{
	int i, w, h;
	int wid = imginf->width;
	int hei = imginf->height;
	int mem=0;
	
	w = wid<<1;
	h = hei<<1;

	//为gauss平滑和dog图像分配空间
	*gaussimg = (float ****)calloc(octaves, sizeof(float ***));
	*dogimg	 = (float ****)calloc(octaves, sizeof(float ***));
	for (i=0; i<octaves; i++)
	{
		w = w>>1;
		h = h>>1;
		mem += mem3D_alloc_float(&(*gaussimg)[i], scales, h, w);
		mem += mem3D_alloc_float(&(*dogimg)[i], scales-1, h, w);		//dog空间比gauss空间少一层
	}

	return mem;
}

/*
释放已分配内存

octave: 尺度空间组数
scales: 每组层数
gaussimg: gauss平滑图像[octaves][scales][height][width]
dogimg:	  dog图像
*/
void gauss_mem_free(int octaves, int scales, float ****gaussimg, float ****dogimg)
{
	mem4D_free_float(gaussimg, octaves, scales);
	mem4D_free_float(dogimg, octaves, scales-1);
}

/*
将图像双线性插值为2倍大小
双线性插值公式:f(i+u,j+v) = (1-u)(1-v)f(i,j) + (1-u)vf(i,j+1) + u(1-v)f(i+1,j) + uvf(i+1,j+1)

imgData:原始图像数据
imginf: 待填充
width:	原始图像宽
height：原始图像高
*/
int imginfo_fill_dbl(imgpix *imgData, ImgInfo **imginf, int width, int height)
{
	int i, j;
	int i1, j1;
	int w = width*2;
	int h = height*2;
	int fij, fij1, fi1j, fi1j1;
	int nv;
	float u,v;

	//为图像分配空间
	*imginf	= (ImgInfo *)calloc(1, sizeof(ImgInfo));
	if (NULL==((*imginf)->imgdata = (imgpix **)calloc(h, sizeof(imgpix *))))
		return 0;
	if (NULL==((*imginf)->imgdata[0]	= (imgpix *)calloc(h*w, sizeof(imgpix))))
		return 0;
	for (i=1; i<h; i++)
		(*imginf)->imgdata[i] = (*imginf)->imgdata[i-1]+w;

	(*imginf)->width		= w;
	(*imginf)->height		= h;
	//插值原始像素
	for (j=0; j<h; j++)
	{
		for (i=0; i<w; i++)
		{
			u = i/2.0; v = j/2.0;
			i1 = i>>1; j1 = j>>1;
			fij		= j1*width + i1;
			fij1	= fij + width;
			fi1j	= fij + 1;
			fi1j1	= fij1 + 1;
			u = u-i1;
			v = v-j1;
			nv	= (int)((1-u)*(1-v)*imgData[fij] + (1-u)*v*imgData[fij1] + u*(1-v)*imgData[fi1j] + u*v*imgData[fi1j1]);
			nv	= (nv<0) ? 0: ((nv<MAX_IMG_VAL) ? nv : MAX_IMG_VAL);
			(*imginf)->imgdata[j][i] = nv;
		}
	}

	return w*h*sizeof(imgpix);
}

/*
填充图像信息

imgData:原始图像数据
imginf: 待填充
width:	原始图像宽
height：原始图像高
*/
int imginfo_fill(imgpix *imgData, ImgInfo **imginf, int width, int height)
{
	int i,j,k;
	(*imginf)				= (ImgInfo *)calloc(1, sizeof(ImgInfo));
	if (((*imginf)->imgdata = (imgpix **)calloc(height, sizeof(imgpix *)))==NULL)
		return 0;
	if (((*imginf)->imgdata[0]	= (imgpix *)calloc(height*width, sizeof(imgpix)))==NULL)
		return 0;
	for (i=1; i<height; i++)
		(*imginf)->imgdata[i] = (*imginf)->imgdata[i-1] + width;

	(*imginf)->width		= width;
	(*imginf)->height		= height;

	for (j=0; j<height; j++)
	{
		for (i=0; i<width; i++)
		{
			k = j*width + i;
			(*imginf)->imgdata[j][i] = imgData[k];
		}
	}

	return height*width*sizeof(imgpix);
}



/*
释放图像

imginfo：要释放的图像结构
*/
void imginfo_release(ImgInfo *imginf)
{
	if (imginf->imgdata)
	{
		if (imginf->imgdata[0])
			free(imginf->imgdata[0]);
		free(imginf->imgdata);
	}
	free(imginf);
}


#ifdef _CRT_SECURE_NO_DEPRECATE
#undef _CRT_SECURE_NO_DEPRECATE
#endif
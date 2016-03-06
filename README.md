#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WIDTH                        1920                  
#define HEIGHT                       1080
#define FRAME                        10                 //Number of sequence
#define v_max                        255
#define BT2020                       1
#define LDR_coeff                    0.6

#define Simple_tonemapping           0

#define Mai_tonemapping              1

#if Mai_tonemapping
#define N                            6
#define Mai_t                        0.3
#define s_max                        231.4
#endif

#define Adaptive_tonemapping         0

#if Adaptive_tonemapping
#define N                            5
#define Adaptive_t                   0.1
#define Max_limit                    1
#define MAX_SLOPE                    300
#define MIN_SLOPE                    150
#endif

typedef unsigned char bitdepth8;
typedef short bitdepth10;

void ChangeYUV444(bitdepth10 *img_in, bitdepth10 *Y_data, bitdepth10 *U_data4, bitdepth10 *V_data4, bitdepth10 *U_data1, bitdepth10 *V_data1, int width, int height, int frame);

void YUV2RGB(bitdepth10 *R_data, bitdepth10 *G_data, bitdepth10 *B_data, bitdepth10 *Y_data, bitdepth10 *U_data4, bitdepth10 *V_data4, int width, int height);

void RGB2Lumin(float *Lumi_data, bitdepth10 *R_data, bitdepth10 *G_data, bitdepth10 *B_data, int width, int height);	

void Lumin2Luma(float *Luma_data, float *Luma_hist, float *Lumi_data, int width, int height);

void getHistogram(int *hist, float *Luma_hist, int width, int height);

#if Mai_tonemapping
void Mai_toneMapping(float *LDR_data, int *hist, float *Luma_data, int width, int height, int n);
#endif

#if Adaptive_tonemapping
void Adaptive_Tonemapping(float *LDR_data, int *hist, float *Luma_data, int width, int height, int n);
#endif

void LDR2RGB(bitdepth8 *LDR_R, bitdepth8 *LDR_G, bitdepth8 *LDR_B, float *LDR_data, float *Lumi_data, bitdepth10 *R_data, bitdepth10 *G_data, bitdepth10 *B_data, int width, int height);

void RGB2YUV(bitdepth8 *LDR_Y, bitdepth8 *LDR_U, bitdepth8 *LDR_V, bitdepth8 *LDR_R, bitdepth8 *LDR_G, bitdepth8 *LDR_B, int width, int height);

int main(void)
{
	FILE *fp_in, *fp_out, *fp_hist;
	bitdepth10 *img_in=0;
	bitdepth8 *img_out=0;
	bitdepth10 *Y_data=0, *U_data1=0, *V_data1=0;   //YUV 4:2:0
	bitdepth10 *U_data4=0, *V_data4=0;              //YUV 4:4:4
	bitdepth10 *R_data=0, *G_data=0, *B_data=0;     //RGB
	float *Lumi_data=0;                             //Luminance
	float *Luma_data=0;                             //Luma
	float *Luma_hist=0;                             //Luma for histogram
	float *LDR_data;                                //LDR Luminance
	bitdepth8 *LDR_R=0, *LDR_G=0, *LDR_B=0;         //LDR RGB
	bitdepth8 *LDR_Y=0, *LDR_U=0, *LDR_V=0;         //LDR YUV

	int x=0, y=0, i=0;
	int num=0;

	fp_in = fopen("Tibul_10.yuv", "rb");
	fp_out = fopen("test_r1.yuv", "wb");

	if(fp_in==NULL)
	{
		printf("File open failed\n");
	}

	//img_in, img_out 메모리 할당
	img_in = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT*1.5*FRAME);
	img_out = (bitdepth8*)malloc(sizeof(bitdepth8)*WIDTH*HEIGHT*1.5*FRAME);

	//YUV 4:2:0 메모리 할당
	Y_data = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);
	U_data1 = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);
	V_data1 = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);

	//YUV 4:4:4 메모리 할당
	U_data4 = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);
	V_data4 = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);

	//RGB 메모리 할당
	R_data = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);
	G_data = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);
	B_data = (bitdepth10*)malloc(sizeof(bitdepth10)*WIDTH*HEIGHT);

	//Luminance, Luma 메모리 할당
	Lumi_data = (float*)malloc(sizeof(float)*WIDTH*HEIGHT);
	Luma_data = (float*)malloc(sizeof(float)*WIDTH*HEIGHT);
	Luma_hist = (float*)malloc(sizeof(float)*WIDTH*HEIGHT);

	//LDR_data 메모리 할당
	LDR_data = (float *)malloc(sizeof(float)*WIDTH*HEIGHT);

	//LDR RGB 메모리 할당
	LDR_R = (bitdepth8 *)malloc(sizeof(bitdepth8)*WIDTH*HEIGHT);
	LDR_G = (bitdepth8 *)malloc(sizeof(bitdepth8)*WIDTH*HEIGHT);
	LDR_B = (bitdepth8 *)malloc(sizeof(bitdepth8)*WIDTH*HEIGHT);

	//LDR YUV 메모리 할당
	LDR_Y = (bitdepth8 *)malloc(sizeof(bitdepth8)*WIDTH*HEIGHT);
	LDR_U = (bitdepth8 *)malloc(sizeof(bitdepth8)*WIDTH*HEIGHT);
	LDR_V = (bitdepth8 *)malloc(sizeof(bitdepth8)*WIDTH*HEIGHT); 

	fread(img_in, sizeof(bitdepth10), WIDTH*HEIGHT*1.5*FRAME, fp_in);
	
	for(num=0; num<FRAME; num++)
	{
		int hist[31] = {0,};
		
		ChangeYUV444(img_in, Y_data, U_data4, V_data4, U_data1, V_data1, WIDTH, HEIGHT, num);

		YUV2RGB(R_data, G_data, B_data, Y_data, U_data4, V_data4, WIDTH, HEIGHT);

		RGB2Lumin(Lumi_data, R_data, G_data, B_data, WIDTH, HEIGHT);

		Lumin2Luma(Luma_data, Luma_hist, Lumi_data, WIDTH, HEIGHT);

		getHistogram(hist, Luma_hist, WIDTH, HEIGHT);

#if 0
		if(num==0)
		{
			fp_hist = fopen("Histogrma.xls", "wt");
			for(i=0; i<31; i++)
			{
				fprintf(fp_hist, "%d\n", hist[i]);
			}
		}
#endif

#if Mai_tonemapping
		Mai_toneMapping(LDR_data, hist, Luma_data, WIDTH, HEIGHT, N);
#endif

#if Adaptive_tonemapping
		Adaptive_Tonemapping(LDR_data, hist, Luma_data, WIDTH, HEIGHT, N);
#endif

		LDR2RGB(LDR_R, LDR_G, LDR_B, LDR_data, Lumi_data, R_data, G_data, B_data, WIDTH, HEIGHT);

		RGB2YUV(LDR_Y, LDR_U, LDR_V, LDR_R, LDR_G, LDR_B, WIDTH, HEIGHT);

		//출력
		for(i=0; i<WIDTH*HEIGHT; i++)
		{
			img_out[(int)(WIDTH*HEIGHT*1.5*num) + i] = LDR_Y[i];
		}
		for(i=0; i<WIDTH*HEIGHT/4; i++)
		{
			img_out[(int)(WIDTH*HEIGHT*1.5*num) + i+(WIDTH*HEIGHT)] = LDR_U[i];
		}
		for(i=0; i<WIDTH*HEIGHT/4; i++)
		{
			img_out[(int)(WIDTH*HEIGHT*1.5*num) + i + (WIDTH*HEIGHT + WIDTH*HEIGHT/4)] = LDR_V[i];
		}
	}

	fwrite(img_out, sizeof(bitdepth8), WIDTH*HEIGHT*1.5*FRAME, fp_out);
	if(fp_out == NULL)
	{
		printf("File open failed\n");
	}

	free(img_in); free(img_out);
	free(Y_data); free(U_data1); free(V_data1);
	free(U_data4); free(V_data4);
	free(R_data); free(G_data); free(B_data);
	free(Lumi_data); free(Luma_data);
	free(LDR_data);
	free(LDR_R); free(LDR_G); free(LDR_B);
	free(LDR_Y); free(LDR_U); free(LDR_V);
}

void ChangeYUV444(bitdepth10 *img_in, bitdepth10 *Y_data, bitdepth10 *U_data4, bitdepth10 *V_data4, bitdepth10 *U_data1, bitdepth10 *V_data1, int width, int height, int frame)
{
	int i=0, j=0;

	for(i=0; i<width*height; i++)
	{
		Y_data[i]=img_in[(int)(width*height*1.5*frame) + i];
	}

	for(i=0; i<(width*height/4); i++)
	{
		U_data1[i]=img_in[(int)(width*height*1.5*frame) + width*height + i];
	}

	for(i=0; i<width*height/4; i++)
	{
		V_data1[i]=img_in[(int)(width*height*1.5*frame) + (width*height)+(width*height/4) +i];
	}

	for(j=0; j<height/2; j++)
	{
		for(i=0; i<width/2; i++)
		{
			U_data4[ (2 * j    ) * width + 2 * i    ] =
			U_data4[ (2 * j    ) * width + 2 * i + 1] =
			U_data4[ (2 * j + 1) * width + 2 * i    ] =
			U_data4[ (2 * j + 1) * width + 2 * i + 1] = U_data1[ j * (width/2) + i ];

			V_data4[ (2 * j    ) * width + 2 * i    ] =
			V_data4[ (2 * j    ) * width + 2 * i + 1] =
			V_data4[ (2 * j + 1) * width + 2 * i    ] =
			V_data4[ (2 * j + 1) * width + 2 * i + 1] = V_data1[ j * (width/2) + i ];
		}
	}
}

void YUV2RGB(bitdepth10 *R_data, bitdepth10 *G_data, bitdepth10 *B_data, bitdepth10 *Y_data, bitdepth10 *U_data4, bitdepth10 *V_data4, int width, int height)
{
	int i=0;
	float Y_temp=0, Cb_temp=0, Cr_temp=0;
	float R_temp=0, G_temp=0, B_temp=0;

	for(i=0; i<width*height; i++)
	{
#if BT2020
		Y_temp = ((float)Y_data[i]-0.0) / 1023.0;
		Cb_temp = ((float)U_data4[i]-512.0) / 1023.0;
		Cr_temp = ((float)V_data4[i]-512.0) / 1023.0;

		R_temp = ((Y_temp) + 1.4746*(float)(Cr_temp))*1023.0;
		G_temp = ((Y_temp) - 0.16455*(float)(Cb_temp) - 0.4592*(float)(Cr_temp))*1023.0;
		B_temp = ((Y_temp) + 1.8814*(float)(Cb_temp))*1023.0;
#else
		R_data[i]= 1.1678*((float)Y_data[i]-(float)64) + 1.6007*((float)V_data4[i]-(float)512);
		G_data[i]= 1.1678*((float)Y_data[i]-(float)64) - 0.3929*((float)U_data4[i]-(float)512) - 0.81532*((float)V_data4[i]-(float)512);
		B_data[i]= 1.1678*((float)Y_data[i]-(float)64) + 2.0232*((float)U_data4[i]-(float)512);
#endif

		R_data[i] = R_temp + 0.5 < 0 ? 0 : R_temp + 0.5;
		R_data[i] = R_temp + 0.5 > 1023 ? 1023 : R_temp + 0.5;
		G_data[i] = G_temp + 0.5 < 0 ? 0 : G_temp + 0.5;
		G_data[i] = G_temp + 0.5 > 1023 ? 1023 : G_temp + 0.5;
		B_data[i] = B_temp + 0.5 < 0 ? 0 : B_temp + 0.5;
		B_data[i] = B_temp + 0.5 > 1023 ? 1023 : B_temp + 0.5;
	}
}

void RGB2Lumin(float *Lumi_data, bitdepth10 *R_data, bitdepth10 *G_data, bitdepth10 *B_data, int width, int height)
{
	int i;

	for(i=0; i<width*height; i++)
	{
#if BT2020
		Lumi_data[i] = 0.2627*(float)(R_data[i])+ 0.6780*(float)(G_data[i]) + 0.0593*(float)(B_data[i]);
#else
		Lumi_data[i] = 0.21258*(float)(R_data[i])+ 0.7154*(float)(G_data[i]) + 0.0721*(float)(B_data[i]);
#endif

		Lumi_data[i] = Lumi_data[i] < 0 ? 0 : Lumi_data[i];
		Lumi_data[i] = Lumi_data[i] > 1023 ? 1023 : Lumi_data[i];
	}
}

void Lumin2Luma(float *Luma_data, float *Luma_hist, float *Lumi_data, int width, int height)
{
	int i;
	float Luma_temp=0;
	int temp=0;

	for(i=0; i<width*height; i++)
	{
		Luma_data[i] = log10((float)Lumi_data[i]);

		temp=Luma_data[i]*10 + 0.5;
		Luma_hist[i] = (float)temp/(float)10;
	}
}

void getHistogram(int *hist, float *Luma_hist, int width, int height)
{
	int i=0, j=0;

	for(i=0; i<width*height; i++)
	{
		j = (int)(10*Luma_hist[i]+0.5);
		hist[j]++;
	}
}

#if Mai_tonemapping
void Mai_toneMapping(float *LDR_data, int *hist, float *Luma_data, int width, int height, int n)
{
	int i=0, j=0;
	int section_num=0, res = 0;
	int *sum=0;
	float *slope=0, *hist_p=0, *y_xis=0;
	float LDR_temp=0;
	float hist_sum=0;
	int max_flag=0;

	sum = (int *)malloc(sizeof(int)*n);
	slope = (float*)malloc(sizeof(float)*n);
	hist_p = (float*)malloc(sizeof(float)*n);
	y_xis = (float*)malloc(sizeof(float)*n);

	section_num = 31/n;
	res = 31%n;

	for(i=0; i<section_num; i++)
	{
		sum[i]=0;
	}
	
	hist[1]+=hist[0];

	for(i=0; i<section_num; i++)
	{
		for(j=1; j<n+1; j++)
		{
			sum[i] += hist[j+(i*n)];
		}
	}

	for(i=0; i<section_num; i++)
	{
		hist_p[i] = (float)sum[i]/(float)(width*height);
	}
	
	for(i=0; i<section_num; i++)
	{
		hist_sum += pow((double)hist_p[i], Mai_t);
	}

	for(i=0; i<section_num; i++)
	{
		slope[i] = ((float)v_max * pow((double)hist_p[i], Mai_t))/((float)n * 0.1 * hist_sum);

		if(slope[i] > s_max)
		{
			max_flag++;
		}

		if(slope[i] <= 0)
		{
			slope[i]= 0; //0.5/(0.1*n);
		}
	}

	if(max_flag > 0)
	{
		for(i=0; i<section_num; i++)
		{
			slope[i] = ((float)v_max * pow((double)hist_p[i], Mai_t))/((float)n * 0.1 * hist_sum);
			
			if(slope[i] > s_max)
			{
				slope[i] = s_max;
			}
			else
			{
				slope[i] = (v_max - (s_max * (float)n * 0.1 * max_flag)) * pow((double)hist_p[i], Mai_t) / ((float)n * 0.1 * hist_sum);
			}
		}
	}

	y_xis[0]=0;
	for(i=1; i<section_num+1; i++)
	{
		y_xis[i] = slope[i-1]*((0.1*n*i) - 0.1*n*(i-1)) + y_xis[i-1];
	}

	for(i=0; i<width*height; i++)
	{
		for(j=0; j<n; j++)
		{
			if((int)(Luma_data[i]*10) >= (int)(0.1*n*j*10 + 1) && (int)(Luma_data[i]*10)<(int)(0.1*n*(j+1)*10 + 1))
			{

				LDR_data[i] = slope[j]*(Luma_data[i] - 0.1*n*j) + y_xis[j];
			}
		}
	}

	for(i=0; i<width*height; i++)
	{
		if((int)(Luma_data[i]*10)>(int)(0) && (int)(Luma_data[i]*10)<(int)(1))
		{
			LDR_data[i] = slope[j]*(Luma_data[i] - 0.1*n*j) + y_xis[j];
		}
	}

#if Simple_tonemapping
	for(i=0; i<width*height; i++)
	{
		LDR_data[i] = (255/3)*Luma_data[i];
	}
#endif
}
#endif

#if Adaptive_tonemapping
void Adaptive_Tonemapping(float *LDR_data, int *hist, float *Luma_data, int width, int height, int n)
{
	int i=0, j=0;
	int *sum=0;
	float *slope=0, *hist_p=0, *y_xis=0, *y_temp, *x_section;
	float *y_dif=0;
	float y_max1=0, y_max2=0, y_max3=0;
	int max_flag=0, min_flag=0;
	int min_num=0;
	float hist_sum=0;

	int hist_n=0;
	int Boundary[N+1];
	int hist_temp=0, temp=0;
	int num=0;

	sum = (int *)malloc(sizeof(int)*n);
	slope = (float*)malloc(sizeof(float)*n);
	hist_p = (float*)malloc(sizeof(float)*n);
	y_xis = (float*)malloc(sizeof(float)*n);
	y_temp = (float*)malloc(sizeof(float)*n);
	x_section = (float*)malloc(sizeof(float)*n);
	y_dif = (float*)malloc(sizeof(float)*n);
	
	hist_n = width * height / n;

	for(i=0; i<31; i++)
	{
		if(hist[i] != 0 && num ==0)
		{
			Boundary[0] = i-1;
			num++;
		}
	}

	num=0;

	for(i=0; i<n; i++)
	{
		sum[i]=0;
	}

	for(i=0; i<31; i++)
	{
		hist_temp += hist[i];

		for(j=0; j<n; j++)
		{
			if(hist_temp >= hist_n * (j+1) && hist_temp < hist_n * (j+2) && num == j)
			{
				Boundary[j+1] = i;
				sum[j] = hist_temp -temp;
				temp += sum[j];
				num++;
			}
		}	
	}

	for(i=0; i<n; i++)
	{
		x_section[i] = Boundary[i+1] - Boundary[i];
	}

	for(i=0; i<n; i++)
	{
		hist_p[i] = (float)sum[i]/(float)(width*height);
	}

	for(i=0; i<n; i++)
	{
		hist_sum += pow((double)hist_p[i], Adaptive_t);
	}

	for(i=0; i<n; i++)
	{
		slope[i] = ((float)v_max * pow((double)hist_p[i], Adaptive_t))/((float)((Boundary[i+1]/10.0) - (Boundary[i]/10.0)) * hist_sum);
	}

#if Max_limit
	y_temp[0]=0;
	for(i=1; i<n+1; i++)
	{
		y_temp[i] = slope[i-1]*((Boundary[i]/10.0) - (Boundary[i-1]/10.0)) + y_temp[i-1];
	}

	for(i=0; i<n; i++)
	{
		if(slope[i] > MAX_SLOPE && y_max1 == 0)
		{
			y_max1 = y_temp[i+1] - (MAX_SLOPE*(x_section[i]/10.0) + y_temp[i]);
			slope[i] = MAX_SLOPE;
			max_flag = 1;
		}
		else if(slope[i] > MAX_SLOPE && max_flag == 1)
		{
			y_max2 = y_temp[i+1] - (MAX_SLOPE*(x_section[i]/10.0) + y_temp[i]);
			slope[i] = MAX_SLOPE;
			max_flag = 2;
		}
		else if(slope[i] > MAX_SLOPE && max_flag == 2)
		{
			y_max3 = y_temp[i+1] - (MAX_SLOPE*(x_section[i]/10.0) + y_temp[i]);
			slope[i] = MAX_SLOPE;
			max_flag = 3;
		}
	}

	for(i=0; i<n; i++)
	{
		if(slope[i] < MIN_SLOPE)
		{
			min_num++;
			min_flag = 1;
		}
	}
	//max_flag == 1 && min_flag == 1
	if(max_flag == 1 && min_flag == 1)
	{
		for(i=0; i<n; i++)
		{
			if(slope[i] < MIN_SLOPE)
			{
				y_dif[i] = y_temp[i+1] - y_temp[i] + y_max1/min_num;
				slope[i] = y_dif[i] / (x_section[i]/10.0);
			}
		}

		y_xis[0]=0;
		for(i=1; i<n+1; i++)
		{
			y_xis[i] = slope[i-1]*((Boundary[i]/10.0) - (Boundary[i-1]/10.0)) + y_xis[i-1];
		}
	}
	//max_flag == 2 && min_flag == 1
	if(max_flag == 2 && min_flag == 1)
	{
		for(i=0; i<n; i++)
		{
			if(slope[i] < MIN_SLOPE)
			{
				y_dif[i] = y_temp[i+1] - y_temp[i] + y_max1/min_num + y_max2/min_num;
				slope[i] = y_dif[i] / (x_section[i]/10.0);
			}
		}

		y_xis[0]=0;
		for(i=1; i<n+1; i++)
		{
			y_xis[i] = slope[i-1]*((Boundary[i]/10.0) - (Boundary[i-1]/10.0)) + y_xis[i-1];
		}
	}
	//max_flag == 3 && min_flag == 1
	if(max_flag == 3 && min_flag == 1)
	{
		for(i=0; i<n; i++)
		{
			if(slope[i] < MIN_SLOPE)
			{
				y_dif[i] = y_temp[i+1] - y_temp[i] + y_max1/min_num + y_max2/min_num + y_max3/min_num;
				slope[i] = y_dif[i] / (x_section[i]/10.0);
			}
		}

		y_xis[0]=0;
		for(i=1; i<n+1; i++)
		{
			y_xis[i] = slope[i-1]*((Boundary[i]/10.0) - (Boundary[i-1]/10.0)) + y_xis[i-1];
		}
	}
#else
	y_xis[0]=0;
	for(i=1; i<n+1; i++)
	{
		y_xis[i] = slope[i-1]*((Boundary[i]/10.0) - (Boundary[i-1]/10.0)) + y_xis[i-1];
	}
#endif
	num=0;
#if 1
	for(i=0; i<width*height; i++)
	{
		for(j=0; j<n; j++)
		{
			if(Luma_data[i] >= Boundary[j]*0.1 && Luma_data[i] < Boundary[j+1]*0.1)
			{

				LDR_data[i] = slope[j]*(Luma_data[i] - 0.1*Boundary[j]) + y_xis[j] < 0 ? 0 : slope[j]*(Luma_data[i] - 0.1*Boundary[j]) + y_xis[j];
				LDR_data[i] = slope[j]*(Luma_data[i] - 0.1*Boundary[j]) + y_xis[j] > 255 ? 255 : slope[j]*(Luma_data[i] - 0.1*Boundary[j]) + y_xis[j];
			}
			else if(Luma_data[i] >= Boundary[n]*0.1)
			{
				LDR_data[i] = slope[n-1]*(Luma_data[i] - 0.1*Boundary[n-1]) + y_xis[n-1] < 0 ? 0 : slope[n-1]*(Luma_data[i] - 0.1*Boundary[n-1]) + y_xis[n-1];
				LDR_data[i] = slope[n-1]*(Luma_data[i] - 0.1*Boundary[n-1]) + y_xis[n-1] > 255 ? 255 : slope[n-1]*(Luma_data[i] - 0.1*Boundary[n-1]) + y_xis[n-1];
			}
		}
	}
	printf("%d\n", num);
#else
	for(i=0; i<width*height; i++)
	{
		for(j=0; j<n; j++)
		{
			if(Luma_data[i] >= Boundary[j]*0.1 && Luma_data[i] < Boundary[j+1]*0.1)
			{
				num=j;
				LDR_data[i] = slope[j]*(Luma_data[i] - 0.1*Boundary[j]) + y_xis[j];

				if(num==3)
				{
					LDR_data[i] = 0;
				}
			}
		}
	}
#endif

#if 0
 	for(i=0; i<width*height; i++)
 	{
 		if(0<=Luma_data[i] && Luma_data[i]<=1.9)
 		{
 			LDR_data[i] = slope[0]*(Luma_data[i]-1.2);
 		}
 		else if(1.9<Luma_data[i] && Luma_data[i]<=2.1)
 		{
 			LDR_data[i] = slope[1]*(Luma_data[i]-1.9) + y_xis[1];
 		}
 		else if(2.1<Luma_data[i] && Luma_data[i]<=2.3)
 		{
 			LDR_data[i] = slope[2]*(Luma_data[i]-2.1) + y_xis[2];
 		}
 		else if(2.3<Luma_data[i] && Luma_data[i]<=2.4)
 		{
 			LDR_data[i] = slope[3]*(Luma_data[i]-2.3) + y_xis[3];
 		}
 		else if(2.4<Luma_data[i] && Luma_data[i]<=3.0)
 		{
 			LDR_data[i] = slope[4]*(Luma_data[i]-2.4) + y_xis[4];
 		}
 		else
 		{
 			printf("에러");
 		}
 	}
#endif
// 	free(sum); free(slope); free(hist_p); free(y_xis);
}
#endif
 
void LDR2RGB(bitdepth8 *LDR_R, bitdepth8 *LDR_G, bitdepth8 *LDR_B, float *LDR_data, float *Lumi_data, bitdepth10 *R_data, bitdepth10 *G_data, bitdepth10 *B_data, int width, int height)
{
	int i=0;

	for(i=0; i<width*height; i++)
	{
		LDR_R[i] = pow((double)(R_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5 < 0 ? 0 : pow((double)(R_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5;
		LDR_R[i] = pow((double)(R_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5 > 255 ? 255 : pow((double)(R_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5;
		LDR_G[i] = pow((double)(G_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5 < 0 ? 0 : pow((double)(G_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5;
		LDR_G[i] = pow((double)(G_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5 > 255 ? 255 : pow((double)(G_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5;
		LDR_B[i] = pow((double)(B_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5 < 0 ? 0 : pow((double)(B_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5;
		LDR_B[i] = pow((double)(B_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5 > 255 ? 255 : pow((double)(B_data[i]/Lumi_data[i]), LDR_coeff) * LDR_data[i] + 0.5;
	}
}

void RGB2YUV(bitdepth8 *LDR_Y, bitdepth8 *LDR_U, bitdepth8 *LDR_V, bitdepth8 *LDR_R, bitdepth8 *LDR_G, bitdepth8 *LDR_B, int width, int height)
{
	int i=0, j=0;
	bitdepth8 *U_temp=0, *V_temp=0;
	float Y_temp2=0, U_temp2=0, V_temp2=0;
	float R_temp=0, G_temp=0, B_temp=0;

	U_temp=(bitdepth8*)malloc(sizeof(bitdepth8)*width*height);
	V_temp=(bitdepth8*)malloc(sizeof(bitdepth8)*width*height);

	for(i=0; i<width*height; i++)
	{
		R_temp = (float)LDR_R[i] / 255.0;
		G_temp = (float)LDR_G[i] / 255.0;
		B_temp = (float)LDR_B[i] / 255.0;

#if BT2020
		Y_temp2 = (0.262700 * R_temp) + (0.678000 * G_temp) + (0.059300 * B_temp);
		U_temp2 = (-0.139630 * R_temp) + (-0.360370 * G_temp) + (0.500000 * B_temp);
		V_temp2 = (0.500000 * R_temp) + (-0.459786 * G_temp) + (-0.040214 * B_temp);
#endif
		Y_temp2 = Y_temp2 * 255.000000;
		U_temp2 = U_temp2 * 255.000000 + 128.000000;
		V_temp2 = V_temp2 * 255.000000 + 128.000000;

		LDR_Y[i] = Y_temp2 + 0.5 < 0 ? 0 : Y_temp2 + 0.5;
		LDR_Y[i] = Y_temp2 + 0.5 > 255 ? 25 : Y_temp2 + 0.5;
		U_temp[i] = U_temp2 + 0.5 < 0 ? 0 : U_temp2 + 0.5;
		U_temp[i] = U_temp2 + 0.5 > 255 ? 255 : U_temp2 + 0.5;
		V_temp[i] = V_temp2 + 0.5 < 0 ? 0 : V_temp2 + 0.5;
		V_temp[i] = V_temp2 + 0.5 > 255 ? 255 : V_temp2 + 0.5;
	}

	for(i=0; i<height/2; i++)
	{
		for(j=0; j<width/2; j++)
		{
			LDR_U[(int)(width/2)*i + j] = U_temp[(int)(width)*(i*2) + (j*2)];
			LDR_V[(int)(width/2)*i + j] = V_temp[(int)(width)*(i*2) + (j*2)];
		}
	}

	free(U_temp); free(V_temp);
}

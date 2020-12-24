/*
*程序功能：适用于非对称预应力布置双端张拉，在参数设置区域和区段要素设置区域输入数据，
*自动计算各区段开始和结束张拉力，区段平均张拉力，张拉力平衡点位置，单端伸长值。
*作者：谢鑫宇
*开发日期：2020.11.23
*/

#include <stdio.h>
#include <math.h>

int main(){

    //参数设置区域：
    //钢绞线区段编号有n段
    int n=9;
    //u：预应力筋与孔道壁的摩擦系数
    float u=0.17;
    //k:孔道每米局部偏差对摩擦的影响系数
    float k=0.0015;
    //钢绞线弹性模量（Mpa）
    float Ep=195000;
    //钢绞线束截面积（mm^2）
    float Ap=2363;
    //张拉控制应力（Mpa）
    float Fc=1339.2;
    //锚下钢绞线张拉力（KN）
    float pss=Ap*Fc/1000;

    //区段要素设置区域：
    //输入区段长度（m）
    float len[]={11.6725,0.7643,21.341,2.042,4.5612,2.042,12.7843,1.5732,1.1482};
    //输入区段曲线段弧度（rad），直线段输入0
    float angle1[]={0,0.095537,0,0.25525,0,0.25525,0,0.19665,0};

    //数据计算区域：
	//区段结构声明
    struct zone{
        double length;
        double angle;
        double ps,pe,pp;
        float elong;
    };
	
    printf("总计有%d个区段，正算编号：1~%d，反算编号：%d~1\n",n,n,n);
    //创建正算结构数组并赋值,计算各区段张拉力、伸长量
    printf("正算数据：\n");
	printf("区段开始张拉力（KN） 区段结束张拉力（KN） 区段平均张拉力（KN） 区段伸长值（m）\n");
    struct zone z[n];
    for (int i = 0; i < n; i++)
    {
        z[i].length=len[i];
        z[i].angle=angle1[i];
        if (i==0)
        {
            z[i].ps=pss;
            z[i].pe=pss*pow(2.7182,(-(k*z[i].length+u*z[i].angle)));
            printf("%d %f %f",i+1,z[i].ps,z[i].pe);
        }else
        {
            z[i].ps=z[i-1].pe;
            z[i].pe=z[i].ps*pow(2.7182,(-(k*z[i].length+u*z[i].angle)));
            printf("%d %f %f",i+1,z[i].ps,z[i].pe);
        }
        z[i].pp=z[i].ps*(1-pow(2.7182,(-(k*z[i].length+u*z[i].angle))))/(k*z[i].length+u*z[i].angle);
        z[i].elong=z[i].pp*z[i].length/(Ep*Ap)*1000;
        printf(" %f %f\n",z[i].pp,z[i].elong);
    }

	//创建反算结构数组并赋值,计算各区段张拉力、伸长量
    printf("反算数据：\n");
	printf("区段开始张拉力（KN） 区段结束张拉力（KN） 区段平均张拉力（KN） 区段伸长值（m）\n");
    struct zone z1[n];
    for (int j = n-1; j > -1; j--)
    {
        z1[j].length=len[j];
        z1[j].angle=angle1[j];
        if (j==n-1)
        {
            z1[j].ps=pss;
            z1[j].pe=pss*pow(2.7182,(-(k*z1[j].length+u*z1[j].angle)));
            printf("%d %f %f",j+1,z1[j].ps,z1[j].pe);
        }else
        {
            z1[j].ps=z1[j+1].pe;
            z1[j].pe=z1[j].ps*pow(2.7182,(-(k*z1[j].length+u*z1[j].angle)));
            printf("%d %f %f",j+1,z1[j].ps,z1[j].pe);
        }
        z1[j].pp=z1[j].ps*(1-pow(2.7182,(-(k*z1[j].length+u*z1[j].angle))))/(k*z1[j].length+u*z1[j].angle);
        z1[j].elong=z1[j].pp*z1[j].length/(Ep*Ap)*1000;
        printf(" %f %f\n",z1[j].pp,z1[j].elong);
    }

	//求解张拉力平衡点在哪个区段，该点距离该段正算起点多少m
    double x; int n_ph;
    for(int l =0; l < n; l++){
		if((z[l].ps-z1[l].pe)*(z[l].pe-z1[l].ps)<=0.0){
            n_ph=l;
			printf("\n张拉力平衡点在第%d段",l+1);
            double a=log(z[l].ps);
            double b=log(z1[l].ps);
            double c=k*z[l].length+u*z[l].angle;
            x=z[l].length*(a-b+c)/(2*c);
            printf("\n张拉力平衡点距该段正算起点%fm\n",x);
		}
	}

	//求解单端伸长值
	//正算单端伸长量
	float sum1=0,sum2=0;
	for(int l =0; l < n_ph; l++){
        sum1=sum1+z[l].elong;
	}
	float pp1=z[n_ph].ps*(1-pow(2.7182,(-(k*x+u*(x/z[n_ph].length)*z[n_ph].angle))))/(k*x+u*(x/z[n_ph].length)*z[n_ph].angle);
    float elong1=pp1*x/(Ep*Ap)*1000;
    printf("正端伸长值：%fm\n",(sum1+elong1));

    //反算单端伸长量
	for(int l =n_ph+1; l < n; l++){
        sum2=sum2+z1[l].elong;
	}
    float x2=z1[n_ph].length-x;
    float pp2=z1[n_ph].ps*(1-pow(2.7182,(-(k*x2+u*(x2/z1[n_ph].length)*z1[n_ph].angle))))/(k*x2+u*(x2/z1[n_ph].length)*z1[n_ph].angle);
    float elong2=pp2*x2/(Ep*Ap)*1000;
    printf("反端伸长值：%fm\n",(sum2+elong2));

    return 0;
}

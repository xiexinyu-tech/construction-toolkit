/*
*�����ܣ������ڷǶԳ�ԤӦ������˫���������ڲ����������������Ҫ�����������������ݣ�
*�Զ���������ο�ʼ�ͽ���������������ƽ����������������ƽ���λ�ã������쳤ֵ��
*���ߣ�л����
*�������ڣ�2020.11.23
*/

#include <stdio.h>
#include <math.h>

int main(){

    //������������
    //�ֽ������α����n��
    int n=9;
    //u��ԤӦ������׵��ڵ�Ħ��ϵ��
    float u=0.17;
    //k:�׵�ÿ�׾ֲ�ƫ���Ħ����Ӱ��ϵ��
    float k=0.0015;
    //�ֽ��ߵ���ģ����Mpa��
    float Ep=195000;
    //�ֽ������������mm^2��
    float Ap=2363;
    //��������Ӧ����Mpa��
    float Fc=1339.2;
    //ê�¸ֽ�����������KN��
    float pss=Ap*Fc/1000;

    //����Ҫ����������
    //�������γ��ȣ�m��
    float len[]={11.6725,0.7643,21.341,2.042,4.5612,2.042,12.7843,1.5732,1.1482};
    //�����������߶λ��ȣ�rad����ֱ�߶�����0
    float angle1[]={0,0.095537,0,0.25525,0,0.25525,0,0.19665,0};

    //���ݼ�������
	//���νṹ����
    struct zone{
        double length;
        double angle;
        double ps,pe,pp;
        float elong;
    };
	
    printf("�ܼ���%d�����Σ������ţ�1~%d�������ţ�%d~1\n",n,n,n);
    //��������ṹ���鲢��ֵ,������������������쳤��
    printf("�������ݣ�\n");
	printf("���ο�ʼ��������KN�� ���ν�����������KN�� ����ƽ����������KN�� �����쳤ֵ��m��\n");
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

	//��������ṹ���鲢��ֵ,������������������쳤��
    printf("�������ݣ�\n");
	printf("���ο�ʼ��������KN�� ���ν�����������KN�� ����ƽ����������KN�� �����쳤ֵ��m��\n");
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

	//���������ƽ������ĸ����Σ��õ����ö�����������m
    double x; int n_ph;
    for(int l =0; l < n; l++){
		if((z[l].ps-z1[l].pe)*(z[l].pe-z1[l].ps)<=0.0){
            n_ph=l;
			printf("\n������ƽ����ڵ�%d��",l+1);
            double a=log(z[l].ps);
            double b=log(z1[l].ps);
            double c=k*z[l].length+u*z[l].angle;
            x=z[l].length*(a-b+c)/(2*c);
            printf("\n������ƽ����ö��������%fm\n",x);
		}
	}

	//��ⵥ���쳤ֵ
	//���㵥���쳤��
	float sum1=0,sum2=0;
	for(int l =0; l < n_ph; l++){
        sum1=sum1+z[l].elong;
	}
	float pp1=z[n_ph].ps*(1-pow(2.7182,(-(k*x+u*(x/z[n_ph].length)*z[n_ph].angle))))/(k*x+u*(x/z[n_ph].length)*z[n_ph].angle);
    float elong1=pp1*x/(Ep*Ap)*1000;
    printf("�����쳤ֵ��%fm\n",(sum1+elong1));

    //���㵥���쳤��
	for(int l =n_ph+1; l < n; l++){
        sum2=sum2+z1[l].elong;
	}
    float x2=z1[n_ph].length-x;
    float pp2=z1[n_ph].ps*(1-pow(2.7182,(-(k*x2+u*(x2/z1[n_ph].length)*z1[n_ph].angle))))/(k*x2+u*(x2/z1[n_ph].length)*z1[n_ph].angle);
    float elong2=pp2*x2/(Ep*Ap)*1000;
    printf("�����쳤ֵ��%fm\n",(sum2+elong2));

    return 0;
}

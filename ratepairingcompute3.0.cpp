/*
本文件和Fp122ByLGWgmp.h结合使用，用于计算 RatePairing,适用于 扩域不可约多项式为 w^12+2,扭曲线为 y^2=x^3+b/w^6的情形
编译命令如下，使用时将 -I 后面改为 Fp12ByLGWgmp.h 的 dirname，需要预先安装大数计算库 -lgmp和其对应的函数库 -lgmpxx;
g++ -g -O2 ratepairingcompute3.0.cpp -o ratepairingcompute3.0 -lgmpxx -lgmp -I ~/Documents/LGWProgrammingTraining/384bitratepairingcompute20180101

*/
// 2018-01-24 edit by LuoGuiwen

#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <Fp12ByLGWgmp.h>
using namespace std;
mpz_class FROBENIUSMATIX[4][6];//用于 G2中的点及 F p12中元素的FrobeniusRaisedto_q

void gTT(Fp12& l,PointInG2Jacobi& DT,const PointInG2Jacobi& T,const PointInG1& P){//正确性已检验
//DT=2*T l=g_{T,T,P}
Fp2 t0,t1,t2,t3,t4,t5,t6,a0,a1;

sqr(t0,T.Xj);
sqr(t1,T.Yj);
sqr(t2,t1);
add(t4,t1,T.Xj);
sqr(t3,t4);
t3-=t0;
t3-=t2;
add(t4,t3,t3);
t3=t4;

mul(t4,t0,3);
add(t6,T.Xj,t4);
sqr(t5,t4);

add(a0,t3,t3);
sub(DT.Xj,t5,a0);

add(a0,T.Yj,T.Zj);
sqr(DT.Zj,a0);
sqr(a0,T.Zj);
DT.Zj-=t1;
DT.Zj-=a0;

sub(a0,t3,DT.Xj);
mul(DT.Yj,a0,t4);
mul(a0,t2,8);
DT.Yj-=a0;

sqr(a0,T.Zj);
mul(a1,t4,a0);
add(a0,a1,a1);
nega(a1,a0);
mul(t3,a1,P.Xa);

sqr(a0,t6);
mul(a1,t1,4);
t6=a0;
t6-=t0;
t6-=t5;
t6-=a1;

sqr(a0,T.Zj);
mul(a1,DT.Zj,a0);
add(t0,a1,a1);
mul(a0,t0,P.Ya);

//l=a0+t3*w+t6*w^3;
l.b0.a0=a0;
l.b1.a0=t3;
l.b1.a1=t6;
}
void gTQ(Fp12& l,PointInG2Jacobi& TAQ,const PointInG2Jacobi& T,const PointInG2& Q,const PointInG1& P){//正确性已验证
//TAQ=T+Q,l=g_{T,Q}(P);
Fp2 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,a0;

sqr(t1,T.Zj);
mul(t0,Q.Xa,t1);
mul(t2,Q.Ya,T.Zj);
add(t3,t2,t2);
sqr(t4,T.Zj);
mul(t1,t3,t4);
sub(t2,t0,T.Xj);
sqr(t3,t2);
mul(t4,t3,4);
mul(t5,t4,t2);

add(t7,T.Yj,T.Yj);
sub(t6,t1,t7);
mul(t9,t6,Q.Xa);
mul(t7,T.Xj,t4);

sqr(t10,t6);
sub(TAQ.Xj,t10,t5);
add(t10,t7,t7);
TAQ.Xj-=t10;

add(t10,T.Zj,t2);
sqr(TAQ.Zj,t10);
sqr(t10,T.Zj);
TAQ.Zj-=t10;
TAQ.Zj-=t3;

sub(t10,t7,TAQ.Xj);
mul(t8,t10,t6);
mul(t10,T.Yj,t5);
add(t0,t10,t10);
sub(TAQ.Yj,t8,t0);

mul(a0,Q.Ya,TAQ.Zj);
add(t10,a0,a0);
add(a0,t9,t9);
sub(t9,a0,t10);
mul(a0,TAQ.Zj,P.Ya);
add(t10,a0,a0);
nega(a0,t6);
t6=a0;

mul(a0,t6,P.Xa);
add(t1,a0,a0);
//l=t10+t1*w+t9*w^3;
l.b0.a0=t10;
l.b1.a0=t1;
l.b1.a1=t9;
}
void FrobeniusRaisedto_q(PointInG2& temp,const PointInG2& Q){//注意我们需要的是G2中的Frobenius变换，不能直接开搞

temp.Xa.x0=Q.Xa.x0*FROBENIUSMATIX[1][2];
temp.Xa.x0%=q;
temp.Xa.x1=Q.Xa.x1*FROBENIUSMATIX[1][2];
temp.Xa.x1%=q;
if(temp.Xa.x1!=0)temp.Xa.x1=q-temp.Xa.x1;

temp.Ya.x0=Q.Ya.x0*FROBENIUSMATIX[1][3];
temp.Ya.x1=Q.Ya.x1*FROBENIUSMATIX[1][3];
temp.Ya.x1%=q;
if(temp.Ya.x1!=0)temp.Ya.x1=q-temp.Ya.x1;

}

void FrobeniusRaisedto_q_1(Fp12& temp,const Fp12& A){
    Fp2 tt;
    conjugate(temp.b0.a0,A.b0.a0);

    conjugate(tt,A.b1.a0);
    mul(temp.b1.a0,tt,FROBENIUSMATIX[1][1]);

    conjugate(tt,A.b0.a1);
    mul(temp.b0.a1,tt,FROBENIUSMATIX[1][2]);

    conjugate(tt,A.b1.a1);
    mul(temp.b1.a1,tt,FROBENIUSMATIX[1][3]);

    conjugate(tt,A.b0.a2);
    mul(temp.b0.a2,tt,FROBENIUSMATIX[1][4]);

    conjugate(tt,A.b1.a2);
    mul(temp.b1.a2,tt,FROBENIUSMATIX[1][5]);
}
void FrobeniusRaisedto_q_2(Fp12& temp,const Fp12& A){
    temp.b0.a0=A.b0.a0;
    mul(temp.b1.a0,A.b1.a0,FROBENIUSMATIX[2][1]);
    mul(temp.b0.a1,A.b0.a1,FROBENIUSMATIX[2][2]);
    mul(temp.b1.a1,A.b1.a1,FROBENIUSMATIX[2][3]);
    mul(temp.b0.a2,A.b0.a2,FROBENIUSMATIX[2][4]);
    mul(temp.b1.a2,A.b1.a2,FROBENIUSMATIX[2][5]);
}
void FrobeniusRaisedto_q_3(Fp12& temp,const Fp12& A){
    Fp2 tt;
    conjugate(temp.b0.a0,A.b0.a0);

    conjugate(tt,A.b1.a0);
    mul(temp.b1.a0,tt,FROBENIUSMATIX[3][1]);

    conjugate(tt,A.b0.a1);
    mul(temp.b0.a1,tt,FROBENIUSMATIX[3][2]);

    conjugate(tt,A.b1.a1);
    mul(temp.b1.a1,tt,FROBENIUSMATIX[3][3]);

    conjugate(tt,A.b0.a2);
    mul(temp.b0.a2,tt,FROBENIUSMATIX[3][4]);

    conjugate(tt,A.b1.a2);
    mul(temp.b1.a2,tt,FROBENIUSMATIX[3][5]);
}

void mulspecial(Fp12& temp,const Fp12& f,const Fp12& l){
//compute the f*l which l has the form l=a+bw+cw^3.
    Fp2 tt1,tt2,tt3,tt4,tt5,tt6;
    Fp6 t1,t2,t3,t4;

    //a=l.b0.a0;b=l.b1.a0;c=l.b1.a1;
    //t2=(b+cv)*f.b1
    mul(tt1,l.b1.a0,f.b1.a0);
    mul(tt2,l.b1.a1,f.b1.a2);
    mul(tt3,l.b1.a1,f.b1.a0);
    mul(tt4,l.b1.a0,f.b1.a1);
    add(t2.a0,tt1,tt2.mulbyu());
    add(t2.a1,tt3,tt4);

    add(tt5,l.b1.a0,l.b1.a1);
    add(tt6,f.b1.a0,f.b1.a1);
    tt6+=f.b1.a2;
    mul(t2.a2,tt5,tt6);
    t2.a2-=tt1;
    t2.a2-=tt2;
    t2.a2-=tt3;
    t2.a2-=tt4;

    //temp=f*l
    mul(t1,f.b0,l.b0.a0);
    add(temp.b0,t1,t2.mulbyv());
    add(t3.a0,l.b0.a0,l.b1.a0);
    t3.a1=l.b1.a1;
    add(t4,f.b0,f.b1);
    mul(temp.b1,t3,t4);
    temp.b1-=t1;
    temp.b1-=t2;
}
void RatePairingFast(Fp12& f,const PointInG2& Q,const PointInG1& P,const mpz_class& a){//正确性已验证，下一步需要写 finalexp
PointInG2Jacobi T,T1;
PointInG2 Q1,Q2;
Fp12 l,l1;

T.Xj=Q.Xa;
T.Yj=Q.Ya;
T.Zj=Id2;

f=Id12;

//准备采用 NAF形式来做
dec2naf(PTNAF,LENGTHNAF,a);//将 e转化为 NAF形式
for(int i=LENGTHNAF-2;i>=0;i--)
{
    gTT(l,T1,T,P);//warning!!最好不要用相同的T,因为一个输入一个输出，有可能会出错！只是这个函数我们检验过不会出错罢了；
    T=T1;
    sqr(l1,f);
    mulspecial(f,l1,l);

    if(NAF[i]<0)
    {
        //Q1=-Q,未考虑Z分量;
        Q1.Xa=Q.Xa;
        nega(Q1.Ya,Q.Ya);

        gTQ(l,T1,T,Q1,P);
        T=T1;
        mulspecial(l1,f,l);
        f=l1;
    }
    else if(NAF[i]>0)
    {
        gTQ(l,T1,T,Q,P);
        T=T1;
        mulspecial(l1,f,l);
        f=l1;
    }

}

FrobeniusRaisedto_q(Q1,Q);
gTQ(l,T1,T,Q1,P);
mulspecial(l1,f,l);
f=l1;
T=T1;

FrobeniusRaisedto_q(Q2,Q1);
Q1.Xa=Q2.Xa;
nega(Q1.Ya,Q2.Ya);
gTQ(l,T1,T,Q1,P);
mulspecial(l1,f,l);
f=l1;
T=T1;



//---finalexp---
Fp12 mp,mp2,mp3,mx,mx2,mx3,mxp,mx2p,mx3p,mx2p2,y5,TT0,TT1;


conjugate(l,f);
inv(l1,f);
mul(f,l,l1);

FrobeniusRaisedto_q_2(l,f);
mul(l1,l,f);
f=l1;

//mpz_class t=mpz_class("29710560942849126597610328210",10);
FrobeniusRaisedto_q_1(mp,f);
FrobeniusRaisedto_q_2(mp2,f);
FrobeniusRaisedto_q_3(mp3,f);
powcyclotomic(mx,f,t);
powcyclotomic(mx2,mx,t);
powcyclotomic(mx3,mx2,t);
FrobeniusRaisedto_q_1(mxp,mx);
FrobeniusRaisedto_q_1(mx2p,mx2);
FrobeniusRaisedto_q_1(mx3p,mx3);
FrobeniusRaisedto_q_2(mx2p2,mx2);

conjugate(y5,mx2);

mul(l,mx3,mx3p);
conjugate(l1,l);
sqrcyclotomic(TT0,l1);

mul(l,mx,mx2p);
conjugate(l1,l);
mul(l,TT0,l1);
mul(TT0,l,y5);

conjugate(l,mxp);
mul(l1,l,y5);
mul(TT1,l1,TT0);

mul(l,TT0,mx2p2);
sqrcyclotomic(l1,TT1);
mul(TT1,l,l1);

sqrcyclotomic(l1,TT1);
conjugate(l,f);
mul(TT0,l1,l);
mul(mx,mp,mp2);
mul(l,mx,mp3);
mul(TT1,l1,l);
sqrcyclotomic(l,TT0);
mul(f,l,TT1);


}

int main(){
//（1）直接调用 C类型的函数类似如下方式；
//mpz_class a,b=5,c=3,d=9;
// mpz_powm(a.get_mpz_t(),b.get_mpz_t(),c.get_mpz_t(),d.get_mpz_t());
//(2)attetion,不要将long, int 等 和 mpz_class做混合运算，会溢出

//-------初始化FROBENIUSMATIX--------
for(int i=1;i<6;i++){
        mpz_class tt=i*(q-1)/12;
        mpz_powm(FROBENIUSMATIX[1][i].get_mpz_t(),beta.get_mpz_t(),tt.get_mpz_t(),q.get_mpz_t());
}

for(int i=1;i<6;i++){
        mpz_class tt=i*(q*q-1)/12;
        mpz_powm(FROBENIUSMATIX[2][i].get_mpz_t(),beta.get_mpz_t(),tt.get_mpz_t(),q.get_mpz_t());
}

for(int i=1;i<6;i++){
        mpz_class tt=i*(q*q*q-1)/12;
        mpz_powm(FROBENIUSMATIX[3][i].get_mpz_t(),beta.get_mpz_t(),tt.get_mpz_t(),q.get_mpz_t());
}


//---------------------------
mpz_class xp1,yp1,xp20,xp21,yp20,yp21;

xp1=mpz_class("1956744441580460502238309352635051631323503114779573815121763847348953982938611\
5078868812251811756396094239887070392",10);
yp1=mpz_class("2757390840807729375166813962558157228161865197849495512520687856347087432516510\
2248436140239416427799806183476526289",10);
xp20=mpz_class("247089144519780884147630805658931622873657265355727504052993336427323186093\
07400903610125295958725276989491791064669",10);
xp21=mpz_class("648904409386212747004167600252328535902692542690034870005012507900253449181804\
396081099851064567786698611373041937",10);
Fp2 X2 (xp20,xp21);

yp20=mpz_class("154525249935744390416116090261784930389872242831498296807494519685727987380\
92806040191166920594537891361189819306619",10);
yp21=mpz_class("719248472050539538379698788804525736093441987505578754439515884560250728545\
8863034418776836442980487120192862632788",10);
Fp2 Y2 (yp20,yp21);
Fp2 Zj0 (1,0);
Fp12 f,g;
PointInG1 P (xp1,yp1);
PointInG2 Q (X2,Y2);
PointInG2 Q1;
PointInG2Jacobi DT,TAQ,T;


cout<<hex;//十六进制输出


cout<<"t:"<<endl<<t<<endl;
cout<<"Trace tr(t):"<<endl<<6*t*t+1<<endl;
cout<<"Characteristic of Ground Field： q(t)"<<endl<<q<<endl;
cout<<"Order of BN Curve E(F_q): n(t)"<<endl<<N<<endl<<endl;

cout<<"Generator of Group G_1: P_1=(x_{P_1},y_{P_1})"<<endl;
cout<<"x_{P_1}："<<endl;
cout<<xp1<<endl;
cout<<"y_{P_1}："<<endl;
cout<<yp1<<endl<<endl;

cout<<"Generator of Group G_2: P_2=(x_{P_2},y_{P_2})"<<endl;
cout<<"x_{P_2}："<<endl;
cout<<xp20<<endl<<xp21<<"*u"<<endl;
cout<<"y_{P_2}："<<endl;
cout<<yp20<<endl<<yp21<<"*u"<<endl<<endl;

//output ratepairing
clock_t start;
int NUM=1000;
start = clock();
for(int i=0;i<NUM;i++)RatePairingFast(f,Q,P,6*t+2);
cout<<"Rate Pairing e(P_1,P_2):"<<endl;
f.oputinseq();
//for(int i=0;i<LENGTHNAF;i++){cout<<*(PTNAF+i)<<" ";}
cout<<endl;
cout<<"384bit-RatePairingComputation Time is \n"<<(double)(clock()-start)/(NUM*CLOCKS_PER_SEC)<<endl;

return 0;
}





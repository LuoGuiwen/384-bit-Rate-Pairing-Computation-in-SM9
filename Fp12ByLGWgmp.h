/*
Only suitable to the tower extension described in Section 4 of the paper.
1-2-6-12 tower extension of F_q_12；
1-2-4-12 tower extension to implement square in F_q_12
*/
// 2018-01-22 edit by LuoGuiwen

#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;
const mpz_class t=mpz_class("29710560942849126597610328210",10);//384bit参数
//const mpz_class t=mpz_class("4593689212103950336",10);//2^62-2^54+2^44;
//const mpz_class t=mpz_class("6917529027646912906",10);//256bit参数
const mpz_class q=36*t*t*t*t+36*t*t*t+24*t*t+6*t+1;//基域和扩域的特征
const mpz_class N=36*t*t*t*t+36*t*t*t+18*t*t+6*t+1;
const mpz_class beta=-2;

int NAF[5000],LENGTHNAF;//转换naf时使用,(q^12-1)/N长度约为 4608c
int *PTNAF=NAF;

//计算一个mpz_class类型数的NAF函数用法，正确性已验证，NAF中存储的数从左到右数位越来越高；
/*
mpz_class bb=132;
dec2naf(PTNAF,LENGTHNAF,bb);
cout<<length<<endl;
*/

inline void dec2naf(int* rpt,int& length,const mpz_class& op){
    // 指针 rpt不能用 int* & rpt形式，因为我们用完指针之后要归位,使得 *PTNAF=NAF,这样才能复用,
mpz_class E=op,rem;
length=0;//length 首先也要归零，不然会出错！
int ee;

    while(E>0)
    {
        if((E%2)==1)
        {
            rem=(E%4);
            ee=mpz_get_ui(rem.get_mpz_t());
            *(rpt+length)=2-ee;
            E-=*(rpt+length);
        }
        else
        {
            *(rpt+length)=0;
        }
        E/=2;
        length+=1;
    }
}

//扩展的欧几里得算法求模逆 u*a mod q=1,换个说法 u=a^(-1) modq
void invxgcd(mpz_class& x1,const mpz_class& a){
x1=1;
mpz_class u=a,v=q,x2=0,quo,x,r;
while(u!=1){
    quo=v/u;
    r=v%u;
    x=quo*x1;
    x=x2-x;
    v=u;
    u=r;
    x2=x1;
    x1=x;
}
x1%=q;
if(x1<0)x1+=q;
}

//binary algorithm for inversion in Fq,长度很大时，速度比xgcd慢
void invBin(mpz_class& x1,const mpz_class& a){
x1=1;
mpz_class u=a,v=q,x2=0;
while(u!=1&&v!=1){
    while(u%2==0){
        u>>=1;
        if(x1%2==0)x1>>=1;
        else{x1+=q;x1>>=1;}
    }
    while(v%2==0){
        v>>=1;
        if(x2%2==0)x2>>=1;
        else{x2+=q;x2>>=1;}
    }
    if(u>=v){u-=v;x1-=x2;}
    else{v-=u;x2-=x1;}
}
if(u==1)x1=x1%q;
else x1=x2%q;

if(x1<0)x1+=q;
}
//上述两个自己实现的算法都不如gmp的这个快呀，下面就用它吧～
//mpz_invert(rop.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
/*
//Barrett reduction ,速度比较慢，耗时3.28e-6;
const mpz_class base=mpz_class("10000000000000000",16);//base=2^64;
const mpz_class barum1=base*base*base*base*base;
const mpz_class barua1=barum1*base*base;
const mpz_class baru=barum1*barua1/q;
void modBar(mpz_class& r,const mpz_class& z){
mpz_class qbar;
qbar=z>>320;
qbar*=baru;
qbar>>=448;

qbar*=q;
qbar%=barua1;
r=z%barua1;
r-=qbar;
if(r<0)r+=barua1;
while(r>=q)r-=q;
}

*/

class Fp2
{
public:
    mpz_class x0,x1;
    Fp2 (){};//直接初始化为 [0,0],必须要有，不然系统不会为无参数的对象调用 constructor.
    Fp2 (const mpz_class& a,const mpz_class& b){x0=a;x1=b;};//一定要输入非负的a,b。为了速度，我们没有检查 a,b的正负
    //gmp中单次运算不会产生临时变量，<gmp-man-6.1.2 page 81>
    //例如 mul(temp.x0,x1,beta)和 temp.x0=x1*beta是一样的直接赋值，所以可以直接用操作符，看起来更简洁；
    inline Fp2 mulbyu(){
        Fp2 temp;temp.x0=x1+x1;
        temp.x0=-temp.x0;
        while(temp.x0<0){temp.x0+=q;}//至多加2次
        temp.x1=x0;return temp;};
    //注意，这里没有 modq
    //u^2-beta=0;return temp是好的，否则会改变对象的值；
    //据说class内自动inline?
    inline void oput(){cout<<"["<<x0<<" "<<x1<<"]"<<endl;};
};

inline void add(Fp2& temp,const Fp2& A,const Fp2& B){//temp=A+B;
    temp.x0=A.x0+B.x0;
    if(temp.x0>=q){temp.x0-=q;};//不必做%,if比while少一次判断
    temp.x1=A.x1+B.x1;
    if(temp.x1>=q){temp.x1-=q;};
}
inline void sub(Fp2& temp,const Fp2& A,const Fp2& B){
    temp.x0=A.x0-B.x0;
    if(temp.x0<0){temp.x0+=q;};//不必做%
    temp.x1=A.x1-B.x1;
    if(temp.x1<0){temp.x1+=q;};//temp=A-B
}
inline void mul(Fp2& temp,const Fp2& A,const Fp2& B){

temp.x0=A.x0*B.x0;
temp.x1=A.x1*B.x1;
temp.x1<<=1;
temp.x0-=temp.x1;//beta=-2
temp.x0%=q;
if(temp.x0<0){temp.x0+=q;};

temp.x1=A.x1*B.x0;
temp.x1+=A.x0*B.x1;
temp.x1%=q;

    //下面的方法虽然在理论上少做一次乘法，但实际实现时反而不如上面的快，
    //究其原因，一是乘法本身很快以至于加法不能忽略，二是引入临时变量非常耗时。
/*
mpz_class s,t;
s=A.x0*B.x0;
t=A.x1*B.x1;
temp.x1=A.x0+A.x1;
temp.x0=B.x0+B.x1;
temp.x1*=temp.x0;
temp.x1-=s;
temp.x1-=t;
//mpz_mod(temp.x1.get_mpz_t(),temp.x1.get_mpz_t(),q.get_mpz_t());
temp.x1%=q;
//if(temp.x1<0){temp.x1+=q;}//if make it take less time.
t<<=1;//-beta=2
temp.x0=s-t;
//mpz_mod(temp.x0.get_mpz_t(),temp.x0.get_mpz_t(),q.get_mpz_t());
temp.x0%=q;
if(temp.x0<0){temp.x0+=q;}
*/
}
//temp=A*b0；b0一定要输入正整数，我们未检查输入；
inline void mul(Fp2& temp,const Fp2& A,const mpz_class& b0){
temp.x0=A.x0*b0;
temp.x0%=q;
temp.x1=A.x1*b0;
temp.x1%=q;};
inline void mulsmallint(Fp2& temp,const Fp2& A,const mpz_class& b0){//乘以小的正整数b0，b0<=4最好
temp.x0=A.x0*b0;
while(temp.x0>=q){temp.x0-=q;}
temp.x1=A.x1*b0;
while(temp.x1>=q){temp.x1-=q;}//罪魁祸首抓出来啦哈哈哈，temp.x1未约减！！！！
}

inline Fp2 operator +(const Fp2& A,const Fp2& B ){Fp2 temp;add(temp,A,B);return temp;};
inline Fp2 operator -(const Fp2& A,const Fp2& B ){Fp2 temp;sub(temp,A,B);return temp;};
inline Fp2 operator *(const Fp2& A,const Fp2& B ){Fp2 temp;mul(temp,A,B);return temp;};
inline Fp2 operator *(const Fp2& A,const mpz_class& b0){Fp2 temp;mul(temp,A,b0);return temp;};
inline Fp2& operator+=(Fp2& x,const Fp2& A) { add(x, x, A); return x;}//这样避免了引入临时变量，厉害！！！；
inline Fp2& operator-=(Fp2& x,const Fp2& A) { sub(x, x, A); return x;}//参考学习 ZZ.h
inline Fp2& operator*=(Fp2& x,const mpz_class& b0) { mulsmallint(x, x, b0); return x;}
//inline Fp2& operator*=(Fp2& x,const Fp2& A) { mul(x, x, A); return x;}//*=要引入临时变量，不然会出错
//inline Fp2& operator*=(Fp2& x,const mpz_class& b0) { mul(x,x,b0); return x;}
//第一个值为输出
inline void sqr(Fp2& temp,const Fp2& A){

    temp.x0=A.x0*A.x0;
    temp.x1=A.x1*A.x1;
    temp.x1<<=1;
    temp.x0-=temp.x1;
    temp.x0%=q;
    if(temp.x0<0){temp.x0+=q;}
    temp.x1=A.x0*A.x1;
    temp.x1<<=1;
    temp.x1%=q;

    /*
    mpz_class t1;
    t1=A.x1+A.x1;
    temp.x1=A.x0-t1;
    t1=A.x0+A.x1;
    temp.x0=t1*temp.x1;
    t1=A.x0*A.x1;
    temp.x0+=t1;
    temp.x0%=q;
    if(temp.x0<0){temp.x0+=q;}

    temp.x1=t1+t1;
    temp.x1%=q;
    */
}

inline Fp2 sqr(const Fp2& A){
    Fp2 temp;sqr(temp,A);return temp;
}

inline void inv(Fp2& temp,const Fp2& A){
    mpz_class t1;
    temp.x0=A.x0*A.x0;
    t1=A.x1*A.x1;
    temp.x1=t1<<1;//t0=t0*(-beta);
    temp.x0+=temp.x1;
    //invxgcd(t1,temp.x0);
    mpz_invert(t1.get_mpz_t(),temp.x0.get_mpz_t(),q.get_mpz_t());
    temp.x0=A.x0*t1;
    temp.x0%=q;//没有检查正负，一定是正的
    temp.x1=A.x1*t1;
    temp.x1=-temp.x1;
    temp.x1%=q;
    if(temp.x1<0)temp.x1+=q;
}
inline void nega(Fp2& temp,const Fp2& A){//temp=-A
    if(A.x0==0)temp.x0=0;
    else temp.x0=q-A.x0;//没有考虑A.x0=0的情形；
    if(A.x1==0)temp.x1=0;
    else temp.x1=q-A.x1;
}
inline void conjugate(Fp2& temp,const Fp2& A){
    temp.x0=A.x0;
    if(A.x1==0)temp.x1=0;
    else temp.x1=q-A.x1;
}

inline void sqrFp4(Fp2& temp0,Fp2& temp1,const Fp2& A0,const Fp2& A1){//Algorithm9,本质上是Fp4中的平方,正确性已验证
    Fp2 t0,t1;
    sqr(t0,A0);
    sqr(t1,A1);
    add(temp0,A0,A1);
    sqr(temp1,temp0);
    temp1-=t0;
    temp1-=t1;

    temp0=t1.mulbyu();
    temp0+=t0;
}


class Fp6//Algorithm12, v^6-beta=0,u^2-beta=0,u=v^3;
{
public:
    Fp2 a0,a1,a2;
    Fp6 (){};//必须要有，不然系统不会为无参数的对象调用 constructor.
    Fp6 (const Fp2& a,const Fp2& b,const Fp2& c){a0=a;a1=b;a2=c;};
    Fp6 mulbyv();
    void oput(){cout<<"["<<endl;a0.oput();a1.oput();a2.oput();cout<<"]"<<endl;};
};
inline Fp6 Fp6::mulbyv(){
    Fp6 temp;
    temp.a0=a2.mulbyu();
    temp.a1=a0;
    temp.a2=a1;
    return (temp);
}
inline void add(Fp6& temp,const Fp6& A,const Fp6& B){add(temp.a0,A.a0,B.a0);add(temp.a1,A.a1,B.a1);add(temp.a2,A.a2,B.a2);};//temp=A+B;
inline void sub(Fp6& temp,const Fp6& A,const Fp6& B){sub(temp.a0,A.a0,B.a0);sub(temp.a1,A.a1,B.a1);sub(temp.a2,A.a2,B.a2);};//temp=A-B
inline void mul(Fp6& temp,const Fp6& A,const Fp6& B){//Algorithm 13写函数虽然快了一丢丢，但是降低了程序的可读性～，哎
/*
//Schoolbook method
    Fp2 t0,t1,t2;
    mul(t0,A.a0,B.a2);
    mul(t1,A.a1,B.a1);
    mul(t2,A.a2,B.a0);
    add(temp.a0,t0,t1);//借用一下temp.a0当做临时变量
    add(temp.a2,temp.a0,t2);//mulbyu()会引入临时变量，故不会出错；

    mul(t0,A.a0,B.a1);
    mul(t1,A.a1,B.a0);
    mul(t2,A.a2,B.a2);
    add(temp.a0,t0,t1);
    add(temp.a1,temp.a0,t2.mulbyu());

    mul(t0,A.a0,B.a0);
    mul(t1,A.a1,B.a2);
    mul(t2,A.a2,B.a1);
    add(temp.a0,t1,t2);
    add(temp.a0,t0,temp.a0.mulbyu());//mulbyu()会引入临时变量，故不会出错；
*/
//两种方法速度相当
//Karatsuba's method
    Fp2 t0,t1,t2,xx;
    mul(t0,A.a0,B.a0);
    mul(t1,A.a1,B.a1);
    mul(t2,A.a2,B.a2);
    add(xx,A.a1,A.a2);
    add(temp.a1,B.a1,B.a2);
    mul(temp.a0,xx,temp.a1);//temp.a1充当临时变量
    temp.a0-=t1;
    temp.a0-=t2;
    temp.a0=temp.a0.mulbyu();
    temp.a0+=t0;

    add(xx,A.a0,A.a1);
    add(temp.a2,B.a0,B.a1);
    mul(temp.a1,xx,temp.a2);
    temp.a1-=t0;
    temp.a1-=t1;
    temp.a1+=t2.mulbyu();

    add(temp.a2,A.a0,A.a2);
    add(xx,B.a0,B.a2);
    temp.a2=temp.a2*xx;
    temp.a2-=t0;
    temp.a2-=t2;
    temp.a2+=t1;
}
inline void mul(Fp6& temp,const Fp6& A,const mpz_class& b0){//temp=A*b0
mul(temp.a0,A.a0,b0);
mul(temp.a1,A.a1,b0);
mul(temp.a2,A.a2,b0);
}
inline void mul(Fp6& temp,const Fp6& A,const Fp2& b0){//temp=A*b0,Algorithm14
mul(temp.a0,A.a0,b0);
mul(temp.a1,A.a1,b0);
mul(temp.a2,A.a2,b0);
}
inline void sqr(Fp6& temp,const Fp6& A){//Algorithm 16

/*
Fp2 c4,c5;

mul(c4,A.a1,A.a2);
add(c5,c4,c4);//能做加法坚决不做乘法，所以这个东东只适用于beta=-2的情形了。
sqr(c4,A.a0);
add(temp.a0,c4,c5.mulbyu());

mul(c4,A.a0,A.a1);
add(c5,c4,c4);
sqr(c4,A.a2);
add(temp.a1,c4.mulbyu(),c5);

mul(c4,A.a0,A.a2);
add(c5,c4,c4);
sqr(c4,A.a1);
add(temp.a2,c4,c5);
*/
//下面的方法少做乘法

Fp2 c3,c4,c5;
add(temp.a0,A.a0,A.a0);
mul(c4,temp.a0,A.a1);
sqr(c5,A.a2);
add(temp.a1,c5.mulbyu(),c4);
sub(temp.a2,c4,c5);
sqr(c3,A.a0);
sub(temp.a0,A.a0,A.a1);
add(c4,temp.a0,A.a2);
add(temp.a0,A.a1,A.a1);
mul(c5,temp.a0,A.a2);
c4=sqr(c4);//不能调用sqr(c4,c4),会出错，一定要这样引入临时变量；
add(temp.a0,c5.mulbyu(),c3);
temp.a2+=c4;
temp.a2+=c5;
temp.a2-=c3;
}
inline void inv(Fp6& temp,const Fp6& A){//Algorithm 17,正确性已检验
Fp2 t0,t1,t2,t3,t4,t5;
sqr(t0,A.a0);
sqr(t1,A.a1);
sqr(t2,A.a2);
mul(t3,A.a0,A.a1);
mul(t4,A.a0,A.a2);
mul(t5,A.a1,A.a2);

t5=t5.mulbyu();
sub(temp.a0,t0,t5);
t2=t2.mulbyu();
sub(temp.a1,t2,t3);
sub(temp.a2,t1,t4);

mul(t2,A.a0,temp.a0);
mul(t5,A.a2,temp.a1);
t5=t5.mulbyu();
t2+=t5;
mul(t5,A.a1,temp.a2);
t5=t5.mulbyu();
t2+=t5;

inv(t5,t2);

mul(t2,temp.a0,t5);
temp.a0=t2;//利用赋值操作避免引入临时变量；
mul(t2,temp.a1,t5);
temp.a1=t2;
mul(t2,temp.a2,t5);
temp.a2=t2;
}
inline void nega(Fp6& temp,const Fp6& A){//temp=-A
nega(temp.a0,A.a0);
nega(temp.a1,A.a1);
nega(temp.a2,A.a2);
}
inline Fp6 operator +(const Fp6& A,const Fp6& B ){Fp6 temp;add(temp,A,B);return temp;};
inline Fp6 operator -(const Fp6& A,const Fp6& B ){Fp6 temp;sub(temp,A,B);return temp;};
inline Fp6 operator *(const Fp6& A,const Fp6& B ){Fp6 temp;mul(temp,A,B);return temp;};
inline Fp6 operator *(const Fp6& A,const mpz_class& B ){Fp6 temp;mul(temp,A,B);return temp;};
inline Fp6 operator *(const Fp6& A,const Fp2& B ){Fp6 temp;mul(temp,A,B);return temp;};
inline Fp6& operator+=(Fp6& x,const Fp6& A) { add(x, x, A); return x;}//这样避免了引入临时变量，厉害！！！；
inline Fp6& operator-=(Fp6& x,const Fp6& A) { sub(x, x, A); return x;}//参考学习 ZZ.h

//inline Fp6& operator*=(Fp6& x,const Fp6& A) { mul(x, x, A); return x;}
//inline Fp6& operator*=(Fp6& x,const mpz_class& b0) { mul(x,x,b0); return x;}
//inline Fp6& operator*=(Fp6& x,const Fp2& b0) { mul(x,x,b0); return x;}

class Fp12
{
public:
    Fp6 b0,b1;
    Fp12 (){};//必须要有，不然系统不会为无参数的对象调用 constructor.
    Fp12 (const Fp6& a,const Fp6& b){b0=a;b1=b;};
    inline void oput(){cout<<"["<<endl;b0.oput();b1.oput();cout<<"]"<<endl;};
    inline void oputinseq(){
    cout<<b0.a0.x0<<"+"<<endl<<b1.a0.x0<<"*w+"<<endl;
    cout<<b0.a1.x0<<"*w^2+"<<endl<<b1.a1.x0<<"*w^3+"<<endl;
    cout<<b0.a2.x0<<"*w^4+"<<endl<<b1.a2.x0<<"*w^5+"<<endl;
    cout<<b0.a0.x1<<"*w^6+"<<endl<<b1.a0.x1<<"*w^7+"<<endl;
    cout<<b0.a1.x1<<"*w^8+"<<endl<<b1.a1.x1<<"*w^9+"<<endl;
    cout<<b0.a2.x1<<"*w^10+"<<endl<<b1.a2.x1<<"*w^11"<<endl;
    };
};

//定义 Fp12中的单位元
const Fp2 OO (0,0);
const Fp2 Id2 (1,0);
const Fp6 omega0 (Id2,OO,OO);
const Fp6 omega1 (OO,OO,OO);
const Fp12 Id12 (omega0,omega1);
/*
//定义 x^12+2的根 w
const Fp2 OO (0,0);
const Fp2 Id (1,0);
const Fp6 omega0 (OO,OO,OO);
const Fp6 omega1 (Id,OO,OO);
const w (omega0,omega1);

//temp=A*w
inline void mulomega(Fp12& temp,const Fp2& A){temp.b1.a0=A;}//其余自动化零即可
inline void mulomega3(Fp12& temp,const Fp2& A){temp.b1.a1=A;}
inline void convFp2toFp12(Fp12& temp,const Fp2& A){temp.b1.a1=A;}
*/

inline void add(Fp12& temp,const Fp12& A,const Fp12& B){add(temp.b0,A.b0,B.b0);add(temp.b1,A.b1,B.b1);};//temp=A+B;
inline void sub(Fp12& temp,const Fp12& A,const Fp12& B){sub(temp.b0,A.b0,B.b0);sub(temp.b1,A.b1,B.b1);};//temp=A-B
inline void mul(Fp12& temp,const Fp12& A,const Fp12& B){//Algorithm20


    Fp6 t0,t1;

    add(t0,A.b0,A.b1);
    add(t1,B.b0,B.b1);
    mul(temp.b0,t0,t1);
    mul(t0,A.b0,B.b0);
    mul(t1,A.b1,B.b1);
    sub(temp.b1,temp.b0,t0);
    temp.b1-=t1;
    add(temp.b0,t0,t1.mulbyv());

//当扩张次数变高以后，所省的一次乘法就非常宝贵啦！，应按照算法来做
/*
    Fp6 c1;
    mul(c1,A.b1,B.b0);
    mul(temp.b0,A.b0,B.b1);//借用temp.b0当临时变量
    add(temp.b1,c1,temp.b0);

    mul(c1,A.b0,B.b0);
    mul(temp.b0,A.b1,B.b1);
    add(temp.b0,c1,temp.b0.mulbyv());
*/
}
inline void sqr(Fp12& temp,const Fp12& A){//Algorithm22,正确性已验证

    //<Faster squaring in the cyclotomic subgroup of six degree extensions>
   /*
    Fp6 t1;
    sqr(t1,A.b1);
    sqr(temp.b0,A.b0);
    t1=t1.mulbyv();
    temp.b0+=t1;
    mul(t1,A.b0,A.b1);
    add(temp.b1,t1,t1);
    */

    Fp6 c2,c3;
    sub(temp.b0,A.b0,A.b1);
    c2=A.b1;
    c2=c2.mulbyv();
    sub(c3,A.b0,c2);
    mul(c2,A.b0,A.b1);
    mul(temp.b1,temp.b0,c3);
    add(temp.b0,temp.b1,c2);
    add(temp.b1,c2,c2);
    c2=c2.mulbyv();
    temp.b0+=c2;
}
inline void sqrcyclotomic(Fp12& temp,const Fp12& A){//正确性已验证，速度非常满意

sqrFp4(temp.b0.a0,temp.b1.a1,A.b0.a0,A.b1.a1);
sqrFp4(temp.b0.a1,temp.b1.a2,A.b1.a0,A.b0.a2);
sqrFp4(temp.b0.a2,temp.b1.a0,A.b0.a1,A.b1.a2);
temp.b1.a0=temp.b1.a0.mulbyu();

temp.b0.a0*=3;
temp.b0.a0-=A.b0.a0;
temp.b0.a0-=A.b0.a0;

temp.b0.a1*=3;
temp.b0.a1-=A.b0.a1;
temp.b0.a1-=A.b0.a1;

temp.b0.a2*=3;
temp.b0.a2-=A.b0.a2;
temp.b0.a2-=A.b0.a2;

temp.b1.a0*=3;
temp.b1.a0+=A.b1.a0;
temp.b1.a0+=A.b1.a0;

temp.b1.a1*=3;
temp.b1.a1+=A.b1.a1;
temp.b1.a1+=A.b1.a1;

temp.b1.a2*=3;
temp.b1.a2+=A.b1.a2;
temp.b1.a2+=A.b1.a2;

}
inline void nega(Fp12& temp,const Fp12& A){//temp is the conjugate of A
nega(temp.b0,A.b0);
nega(temp.b1,A.b1);
}

inline void conjugate(Fp12& temp,const Fp12& A){//temp is the conjugate of A
temp.b0=A.b0;
nega(temp.b1,A.b1);
}

inline void inv(Fp12& temp,const Fp12& A){//Algorithm23,正确性已检验
    //速度还可以优化，但鉴于在求pairing的过程中只需要算一次逆元，可以将就用。而且，速度已经比ntl快啦
Fp6 t0,t1;
sqr(t0,A.b0);
sqr(t1,A.b1);
t1=t1.mulbyv();
t0-=t1;
inv(t1,t0);
mul(temp.b0,A.b0,t1);
nega(t0,A.b1);
mul(temp.b1,t0,t1);
}

//仅用于验证算法正确性，不在程序中使用
void pow(Fp12& temp,const Fp12& A,const mpz_class& e){//temp=A^e,正确性已检验,速度还可以，比预想的快
temp=A;
Fp12 temp1,invA;
inv(invA,A);
dec2naf(PTNAF,LENGTHNAF,e);
for(int i=LENGTHNAF-2;i>=0;i--)
{
    sqr(temp1,temp);
    temp=temp1;
    if(NAF[i]!=0)
        {
            if(NAF[i]>0)
                {
                    mul(temp1,temp,A);
                    temp=temp1;
                }
            else
            {
                mul(temp1,temp,invA);
                temp=temp1;
            }
        }
}
}

void powcyclotomic(Fp12& temp,const Fp12& A,const mpz_class& e){//temp=A^e,Algorithm 25,

temp=A;
Fp12 temp1,invA;//引入两个临时变元
conjugate(invA,A);//invA=A^(-1)
dec2naf(PTNAF,LENGTHNAF,e);//将 e转化为 NAF形式
for(int i=LENGTHNAF-2;i>=0;i--)
{
    sqrcyclotomic(temp1,temp);
    temp=temp1;
    if(NAF[i]!=0)
        {
            if(NAF[i]>0)
                {
                    mul(temp1,temp,A);
                    temp=temp1;
                }
            else
            {
                mul(temp1,temp,invA);
                temp=temp1;
            }
        }
}

}

class PointInG1
{
public:
    mpz_class Xa,Ya;//仿射坐标
    PointInG1 (){};//必须要有，不然系统不会为无参数的对象调用 constructor.
    PointInG1 (const mpz_class& a,const mpz_class& b){Xa=a;Ya=b;};//未检查大小，初始化时一定要注意0=<a,b<q
};

class PointInG2
{
public:
    Fp2 Xa,Ya;//仿射坐标
    PointInG2 (){};//必须要有，不然系统不会为无参数的对象调用 constructor.
    PointInG2 (const Fp2& a,const Fp2& b){Xa=a;Ya=b;};
};

class PointInG2Jacobi
{
public:
    Fp2 Xj,Yj,Zj;//Jacobian加重射影坐标
    PointInG2Jacobi (){};//必须要有，不然系统不会为无参数的对象调用 constructor.
    PointInG2Jacobi (const Fp2& a,const Fp2& b,const Fp2& c){Xj=a;Yj=b;Zj=c;};//未检查大小，初始化时一定要注意0=<a,b<q
};

//-------------------------




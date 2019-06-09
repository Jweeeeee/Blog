# 对伊藤微分公式和Black-Scholes公式的理解


## 0.引言 

经过两周对布朗运动这一章的学习(由于中间穿插了第六章，所以整体时间更长)，也尝试过常微分方程和常系数随机微分方程的计算机模拟，又阅读了几篇文章，对其中的重点：伊藤微分公式和一般随机微分方程有一些理解，在这里总结一下。

## 1.布朗运动

### 1.1 定义

定义：随机过程${X(t),t\geq0}$称为Brown运动，如果它满足如下三个条件：
$$
(i)X(0)=0；\\
(ii)随机过程X有平稳独立增量；\\
(iii)对每个t>0，X（t)服从N(0,C^2t)
$$


若$c=1$则称其为标准Brown运动。

从定义我们可以知道：

1.标准Brown运动在$t=0$时的状态为$0$；

2.可以推出Brown运动是一个**马尔科夫过程**，任意$t$时刻之后的状态仅和$t$时刻的状态有关，而与历史无关，另外还可以证明它是**鞅过程**和**正态过程**（即**高斯过程**）；

3.在任何有限时间区间内标准Brown运动的变化服从均值为0，方差为${\Delta}t$的正态分布$N(0,{\Delta}t)$，**而且其方差会随着时间区间长度线性增加**。

### 1.2 性质

#### 1.1.2.1 不变性

一些简单的不变性列举，其证明大多是利用**平稳独立增量过程**这个性质。

若$B_t$是**标准Brown运动**，则

**(1)对称性：$Y_t=-B_t$**

**(2)起点变换：$Y_t=B_{t+s}-B_s$**

**(3)尺度变换：$Y_t=cB_\frac{t}{c^2}$**

**(4)时间倒置：$Y_t=tB_{\frac{1}{t}}$**

**(5)时间反向：$Y_t=B_u-B_{u-t}$**

也是**标准Brown运动**。

#### 1.1.2.2 Brown处处不可微

对于任给的正数M，有
$$
P(|\frac{X(t_0+{\Delta}t)-X(t_0)}{{\Delta}t}|\leq M)\\
=P（|\frac{X（{\Delta}t}{\sqrt{{\Delta}t}}）|）\\
=\phi(M\sqrt{{\Delta}t})-\phi(-M\sqrt{{\Delta}t})\to 0({\Delta}t \to 0)
$$
这可能是最好理解的性质：Brown运动是连续的，但它在任一点$t_0$的导数有限的概率为0，i.e，对几乎每条样本轨道上任意一点$t_0$，其导数不存在，也就是说**固定$t_0$，Brown运动不可导**。进一步可以证明**Brown运动处处不可微**（证明没啃清白）。

#### 1.1.2.3 其它性质

对书上其他的性质理解不是很深，所以来说一下在别的地方看到的性质。

**(1)**Brown运动的轨迹会频繁的穿越时间轴$t$，即在时间轴上下波动，这一点其实就是书上对**Brown运动每个状态$a$都常返（a是零常返）**的证明

**(2)**在任意时刻$t$，它的位置$B(t)$不会偏离**正负一个标准差（$B(0)\pm\sqrt{\Delta t}$）**太远

## 2.$(dW(t))^2=dt$的导出

### 2.1  连续可微函数$f(t)$的二次变分

这个概念从别的地方看的，书上只讲了**Brown运动的二次变差过程**，也就是$(dW(t))^2=\displaystyle \lim_{\lambda_n \to 0}\sum^{n}_{i=1}(W(t_i)-W(t_{i-1}))^2$

**定义：**

考虑时间区间$[0,T]$和该区间内的一个划分 ，$\Pi=\{0=t_0<t_1<t_2...<t_N=T\}$则对于任意一个连续函数 $f(t)$，它的二次变分（quadratic variation）定义为：
$$
\sum^{N-1}_{i=0}[f(t_{i+1})-f(t_i)]^2
$$
**推论：**

对于一个连续且在$[0,T]$上处处可微的函数$f(t)$，可以由中值定理得出$\sum^{N-1}_{i=0}[f(t_{i+1})-f(t_i)]^2\leq max_{s\in[0,T]}f'(s)^2\cdot max\{t_{i-1}-t_i\}\cdot T$

由此，对区间$[0,T]$分割足够细时，$||\Pi||=max\{t_{i+1}-t_i\}\to 0$，函数$f(t)$的二次变分为$0$

### 2.2 Brown运动的二次变分

把上述$f(t)$换成$B(t)$即可，Brown运动的二次变分：
$$
\sum^{N-1}_{i=0}[B(t_{i+1})-B(t_i)]^2
$$
但推论有变化:
$$
\displaystyle \lim_{||\Pi|| \to 0}\sum^{n}_{i=1}(B(t_i+1)-B(t_{i}))^2=T
$$
即，对区间$[0,T]$分割足够细时，$||\Pi||=max\{t_{i+1}-t_i\}\to 0$，随机过程$B(t)$的二次变分为$T$（区间长度），而不是0

**理解：**

对于**Brown运动**，其非零的二次变分说明**随机性使得它的波动太频繁**，以至于不管我们如何细分区间$T$、得到多么微小的划分区间，这些微小区间上的**位移差的平方逐段累加起来的总和(二次变分的几何意义)**都不会消失（即二次变分不为0），而是等于这个**区间的长度 $T$**

### 2.3  $(dW(t))^2=dt$

综上，Brown运动的二次变分公式也可以写成$(dB)^2=dt$，这是**伊藤微分公式**推导的关键。

### 2.4 $dW(t)=\sqrt{dt}Z,Z服从N(0,1)$

如何理解这个式子呢？先将其写成增量的形式：
$$
{\Delta}W(t)=\sqrt{dt}Z
$$
对比一般的确定性函数$f$增量和微分的关系：
$$
{\Delta}f(t)=f'(t)dt+o(dt)
$$
我们发现Brown运动的增量与$\sqrt{\Delta t}$成正比，与一般的确定性函数$f$增量和微分的关系不同的是，**Brown运动的增量和微分不再具有线性关系**，也就表明在Brown的样本轨道的任意一点附近不能“以直代曲”。这也构成了随机微分方程和确定性微分方程的本质区别。

## 3.多元函数的泰勒展开

若函数$f$在点$P_0(x_0,y_0)$的某领域$U(P_0)$上有直到$n+1$阶的连续偏导数，则对$U(P_0)$内任一点$(x_0+h,y_0+k)$，存在相应的${\theta}\in(0,1)$，使得
$$
f(x_0+h,y0+k) =
f(x_0,y_0) +
(h\frac{\partial }{\partial x}+k\frac{\partial }{\partial y})f(x_0,y_0) + \\
\frac{1}{2!}(h\frac{\partial }{\partial x}+k\frac{\partial }{\partial y})^2f(x_0,y_0)+\cdots + \\
\frac{1}{n!}(h\frac{\partial }{\partial x}+k\frac{\partial }{\partial y})^nf(x_0+{\theta}h,y_0+{\theta}k)
$$


其中，
$$
(h\frac{\partial }{\partial x}+k\frac{\partial }{\partial y})^mf(x_0+{\theta}h,y_0+{\theta}k) =
\sum^{m}_{i=0}C^i_m\frac{\partial }{\partial x^i \partial y^{m-i}}f(x_0,y_0)h^ik^{m-i}\\
h={\Delta}x,k={\Delta}y
$$
若只需求$R_n=o(\rho^n)(\rho=\sqrt{{\Delta}x^2+{\Delta}y^2})$，则只需$f$在$U(P_0)$内存在直到$n$阶连续偏导数，便有
$$
f(x_0+{\Delta}x,y_0+{\Delta}y)\\
=\\
f(x_0,y_0) +
\sum^{n}_{P=1} \frac{1}{p!}({\Delta}x\frac{\partial }{\partial x} +{\Delta}y\frac{\partial }{\partial y})^Pf(x_0,y_0) +
o(\rho^n)
$$
这个公式将帮助我们导出**伊藤微分公式**

## 4.伊藤微分公式

### 4.1 伊藤微分公式

$Ito微分公式$ 设实函数$f(x,y)$关于$x$有二阶连续偏导数，关于$y$有一阶连续偏导数，若${W(t),t\geq 0}$是参数为$\sigma^2$的Brown运动，则
$$
{\Delta}f(W(t),t) =
\frac{\partial f}{\partial x}dW(t) +
(\frac{\partial f}{\partial y} +
\frac{\sigma^2}{2}\frac{\partial^2 f}{\partial x^2})dt +
o(dt)
$$
书上给出的证明条件是$f$关于$x$和$y$都有二阶连续偏导数。

证明思路是对$f(x,y)$进行泰勒展开，展到二阶，然后处理掉其中的无穷小项。具体过程就不摆了，简单的写一下思路以及理解了的点吧。

**(1)从$df=\big(\frac{dBt}{dt}f'(Bt)\big)dt$到$df=f'(Bt)dBt$**

前者显然是直观的微分形式，但由于Brown运动处处不可导，所以这样的微分是不可行的；

后者绕开了$\frac{dBt}{dt}$，但是这样也是错误的，这是由于**Brown运动的二次变分非零**。当我们用泰勒展开写出它的前两项时，就明白为什么后者也是不可行了。

**(2)要展开到二阶的原因**

由一般函数的泰勒展开中：
$$
f(x+\Delta x)-f(x)=f'(x)\Delta x +
\frac{f''(x)}{2!}\Delta x^2 +
\frac{f'''(x)}{3!}\Delta x^3 \cdots +
\frac{f^n(x)}{n!}\Delta x^n
$$
从第二项开始$\Delta x^m（m>1）$都是$\Delta x$的**高阶无穷小**，所以可以略去，只留第一项，

而**Brown运动**则不行，二阶偏导会出现$(dBt)^2=dt$，不再是高阶无穷小，所以**无法略去**；

**(3)无穷小项的处理**

$(dW(t))^2=\sigma^2dt+o(dt)$，$dW(t)dt=\sigma(dt+o(dt)^\frac{3}{2})$，$dt^2=o(dt)$，第三个显然，第一个和第二个用到了前面的**2.3**和**2.4**。

### 4.2 一般随机微分方程

扩散方程模型：
$$
dX(t)=\mu \big(X(t),t\big)dt + \sigma \big(X(t),t\big)dW(t)
$$
其中$\mu \big(X(t),t\big)$和$\sigma \big(X(t),t\big)$是$X(t)$和$t$的函数。

令$Y(t)=f(X(t),t)$，推导随机过程$Y$满足的随机微分方程：
$$
\Delta Y(t) =
\frac{\partial f}{\partial x}dX(t) +
\frac{\partial f}{\partial y}dt +
\frac{1}{2}
\big(
\frac{\partial^2 f}{\partial x^2}(dX(t))^2 +
2\frac{\partial^2f}{\partial x \partial t}dX(t)dt +
\frac{\partial^2f}{\partial t^2}(dt)^2
\big) +
o(dt)
$$
将$dX(t)=\mu \big(X(t),t\big)dt + \sigma \big(X(t),t\big)dW(t)$代入上面方程，其中，
$$
(X(t))^2 = \mu^2dt^2 + \sigma^2(dW(t))^2 + 2\mu\sigma dW(t)dt = \sigma^2dt + o(dt) \\
dXtdt = \mu dt^2 + \sigma dW(t)dt = o(dt) \\
(dt)^2 = o(dt)
$$
忽略高阶无穷小项，可得：
$$
dY(t) =
\big(
\frac{\partial f}{\partial x}\mu(X(t),t) +
\frac{\partial f}{\partial t} + 
\frac{1}{2}\frac{\partial^2 f}{\partial x^2}\sigma^2(X(t),t)
\big)dt +
\frac{\partial f}{\partial x}\sigma(X(t),t)dW(t)
$$
从这里也可以感受到随机微分方程的解往往是先猜解后验证。

### 4.3 几何布朗运动

设随机过程$S$满足
$$
dS(t) = \mu S(t) dt + \sigma S(t) dW(t)
$$
其中$\mu，\sigma>0$为常数，$W(t)$为标准Brown运动，满足上述微分方程的解称为几何Brown运动。

在这里给出其解：
$$
S(t) = exp\{\sigma W(t) + \big(\mu - \frac{\sigma^2}{2}t\big)\}
$$

## 5.Black-Scholes公式

这里省略介绍$BS$公式的经济学背景，从数学上看，$Black-Scholes$公式其实就是在思考如何消除$\Delta W(t)$。

$f=f(S(t),t)$满足SDE：
$$
df = 
df(S(t),t) =
\big(
\frac{\partial f}{\partial x}\mu S +
\frac{\partial f}{\partial t} + 
\frac{1}{2}\frac{\partial^2 f}{\partial x^2}\sigma^2 S^2
\big)dt +
\frac{\partial f}{\partial x}\sigma S dW(t)
$$
$S=S(t)$满足SDE：
$$
dS = \mu S {\Delta} t + \sigma S {\Delta} t
$$
定义证券组合价值为$\Pi$，其满足：
$$
\Pi = -f + \frac{\partial f}{\partial s} \cdot S \\
d\Pi = -df + \frac{\partial f}{\partial s} \cdot dS
$$
将$df = df(S(t),t)$和$dS = \mu S {\Delta} t + \sigma S {\Delta} t$代入上式，可得：
$$
d\Pi = 
\big(
-\frac{\partial f}{\partial t} -
\frac{1}{2}\frac{\partial^2 f}{\partial^2 t}\sigma^2S^2
\big)dt
$$
这里$dW(t)$被抵消掉了，也就是消去了瞬时收益率的风险项。

在不存在无风险套利的市场中，该投资组合的瞬时收益率$d\Pi$必须等于无风险收益率$r$，即
$$
d\Pi = r \Pi dt
$$
将$d\Pi = 
\big(
-\frac{\partial f}{\partial t} -
\frac{1}{2}\frac{\partial^2 f}{\partial^2 t}\sigma^2S^2
\big)dt$和$\Pi = -f + \frac{\partial f}{\partial s} \cdot S$代入上式，可得：
$$
-\big(
\frac{\partial f}{\partial t} + 
\frac{1}{2}\frac{\partial^2 f}{\partial s^2}\sigma^2S^2
\big)dt =
r\big( 
\frac{\partial f}{\partial s}S -
f
\big)dt
$$
化简得：
$$
\frac{\partial f}{\partial t} +
r S \frac{\partial f}{\partial s} + 
\frac{1}{2}\frac{\partial^2 f}{\partial s^2}\sigma^2 S^2 =
rf
$$
上式称为$Black-Scholes$微分方程。

**[参考资料]**

[《随机过程 方兆本 第三版》](https://max.book118.com/html/2019/0309/8043000015002012.shtm)

[布朗运动、伊藤引理、BS公式（前篇）](https://zhuanlan.zhihu.com/p/38293827)

[布朗运动、伊藤引理、BS公式（后篇）](https://zhuanlan.zhihu.com/p/38294971)

[经济金融系列学习：伊藤引理](<https://www.jianshu.com/p/d7abb9e1ed7d>)

## 曲线拟合数学模型建立
### 例1
$y=ap^2+bp+c$
#### 定义观测值
$p_i,y_i$
#### 定义待估计变量
$x=\begin{bmatrix}a\\b\\c\end{bmatrix}$
#### 设计误差模型
$\hat{y_i}=e^{ap_i^2+bp_i+c}$

$e_i=y_i-\hat{y_i}=y_i-e^{ap_i^2+bp_i+c}$
#### 求解误差项对待估计变量的偏导数
$\frac{\partial e_i}{\partial x}=\frac{\partial e_i}{\partial \hat{y_i}}\frac{\partial \hat{y_i}}{\partial x}$

$\frac{\partial e_i}{\partial \hat{y_i}}=\frac{\partial (y_i-\hat{y_i})}{\partial \hat{y_i}}=-1$

$\frac{\partial \hat{y_i}}{\partial x}=\begin{bmatrix}\frac{\partial \hat{y_i}}{\partial a}\\\frac{\partial \hat{y_i}}{\partial b}\\\frac{\partial \hat{y_i}}{\partial c}\end{bmatrix}$

$\frac{\partial \hat{y_i}}{\partial a}=\frac{\partial (e^{ap_i^2+bp_i+c})}{\partial a}=p_i^2e^{ap_i^2+bp_i+c}=p_i^2\hat{y_i}$

$\frac{\partial \hat{y_i}}{\partial b}=\frac{\partial (e^{ap_i^2+bp_i+c})}{\partial b}=p_ie^{ap_i^2+bp_i+c}=p_i\hat{y_i}$

$\frac{\partial \hat{y_i}}{\partial c}=\frac{\partial (e^{ap_i^2+bp_i+c})}{\partial c}=e^{ap_i^2+bp_i+c}=\hat{y_i}$

$J_i=\frac{\partial e_i}{\partial x}=\frac{\partial e_i}{\partial \hat{y_i}}\frac{\partial \hat{y_i}}{\partial x}=-\begin{bmatrix}p_i^2\hat{y_i}\\p_i\hat{y_i}\\\hat{y_i}\end{bmatrix}=-\hat{y_i}\begin{bmatrix}p_i^2\\p_i\\1\end{bmatrix}$
### 例2
$y = a^p + p^b - sin(cp) + e^{ap^2 - bp + \frac{c}{p}} - ln(dp+e)$
#### 定义观测值
$p_i,y_i$
#### 定义待估计变量
$x=\begin{bmatrix}a\\b\\c\\d\\e\end{bmatrix}$
#### 设计误差模型
$\hat{y_i}=a^{p_i}+p_i^b-sin(cp_i)+e^{ap_i^2-bp_i+\frac{c}{p_i}}-ln(dp_i+e)$

$e_i=y_i-\hat{y_i}=y_i-a^{p_i}-p_i^b+sin(cp_i)-e^{ap_i^2-bp_i+\frac{c}{p_i}}+ln(dp_i+e)$
#### 求解误差项对待估计变量的偏导数
$\frac{\partial e_i}{\partial x}=\frac{\partial e_i}{\partial \hat{y_i}}\frac{\partial \hat{y_i}}{\partial x}$

$\frac{\partial e_i}{\partial \hat{y_i}}=\frac{\partial (y_i-\hat{y_i})}{\partial \hat{y_i}}=-1$

$\frac{\partial \hat{y_i}}{\partial x}=\begin{bmatrix}\frac{\partial \hat{y_i}}{\partial a}\\\frac{\partial \hat{y_i}}{\partial b}\\\frac{\partial \hat{y_i}}{\partial c}\\\frac{\partial \hat{y_i}}{\partial d}\\\frac{\partial \hat{y_i}}{\partial e}\end{bmatrix}$

$\frac{\partial \hat{y_i}}{\partial a}=\frac{\partial (a^{p_i}+p_i^b-sin(cp_i)+e^{ap_i^2-bp_i+\frac{c}{p_i}}-ln(dp_i+e))}{\partial a}=p_ia^{p_i-1}+p_i^2e^{ap_i^2-bp_i+\frac{c}{p_i}}$

$\frac{\partial \hat{y_i}}{\partial b}=\frac{\partial (a^{p_i}+p_i^b-sin(cp_i)+e^{ap_i^2-bp_i+\frac{c}{p_i}}-ln(dp_i+e))}{\partial b}=p_i^bln(p_i) - p_ie^{ap_i^2 - bp_i + \frac{c}{p_i}}$

$\frac{\partial \hat{y_i}}{\partial c}=\frac{\partial (a^{p_i}+p_i^b-sin(cp_i)+e^{ap_i^2-bp_i+\frac{c}{p_i}}-ln(dp_i+e))}{\partial c}=-p_icos(cp_i) + \frac{1}{p_i}e^{ap_i^2 - bp_i + \frac{c}{p_i}}$

$\frac{\partial \hat{y_i}}{\partial d}=\frac{\partial (a^{p_i}+p_i^b-sin(cp_i)+e^{ap_i^2-bp_i+\frac{c}{p_i}}-ln(dp_i+e))}{\partial d}=-\frac{p_i}{dp_i+e}$

$\frac{\partial \hat{y_i}}{\partial e}=\frac{\partial (a^{p_i}+p_i^b-sin(cp_i)+e^{ap_i^2-bp_i+\frac{c}{p_i}}-ln(dp_i+e))}{\partial e}=-\frac{1}{dp_i+e}$

$J_i=\frac{\partial e_i}{\partial x}=\frac{\partial e_i}{\partial \hat{y_i}}\frac{\partial \hat{y_i}}{\partial x}=-\begin{bmatrix}p_ia^{p_i-1}+p_i^2e^{ap_i^2-bp_i+\frac{c}{p_i}}\\p_i^bln(p_i) - p_ie^{ap_i^2 - bp_i + \frac{c}{p_i}}\\-p_icos(cp_i) + \frac{1}{p_i}e^{ap_i^2 - bp_i + \frac{c}{p_i}}\\-\frac{p_i}{dp_i+e}\\-\frac{1}{dp_i+e}\end{bmatrix}=\begin{bmatrix}-p_ia^{p_i-1}-p_i^2e^{ap_i^2-bp_i+\frac{c}{p_i}}\\-p_i^bln(p_i) + p_ie^{ap_i^2 - bp_i + \frac{c}{p_i}}\\p_icos(cp_i) - \frac{1}{p_i}e^{ap_i^2 - bp_i + \frac{c}{p_i}}\\\frac{p_i}{dp_i+e}\\\frac{1}{dp_i+e}\end{bmatrix}$

## PnP问题
### 定义观测值
三维空间中的点：$P_i=\begin{bmatrix}X_i\\Y_i\\Z_i\end{bmatrix}$

二维图像中的点：$u_i=\begin{bmatrix}u_i\\v_i\end{bmatrix}$
### 待估计变量
位姿：

$\mathfrak{se}(3)=\left\{\xi=\begin{bmatrix}\rho\\\phi\end{bmatrix}\in \mathbb{R}^6, \rho\in \mathbb{R}^3, \phi\in \mathfrak{so}(3),\xi^{\wedge}=\begin{bmatrix}\phi^{\wedge}&\rho\\\mathbf{0}^T&0\end{bmatrix}\in\mathbb{R}^{4\times 4}\right\}$

$SE(3)=\left\{T=\begin{bmatrix}R&t\\0^T&1\end{bmatrix}\in \mathbb{R}^{4\times 4}, R\in SO(3), t\in \mathbb{R}^3\right\}$
### 误差模型
$s_i\begin{bmatrix}u_i\\v_i\\1\end{bmatrix}=KTP_i=\begin{bmatrix}f_x&0&c_x\\0&f_y&c_y\\0&0&1\end{bmatrix}\begin{bmatrix}R&t\\0^T&1\end{bmatrix}\begin{bmatrix}X_i\\Y_i\\Z_i\\1\end{bmatrix}$

$\hat{u_i}=\frac{1}{s_i}KTP_i$

$T^*=\arg \underset{T}{min}\frac{1}{2}\mathop{\sum}\limits_{i=1}^n\left\|u_i-\hat{u_i}\right\|^2$

$e=u-\hat{u}$

相机坐标系下空间点的坐标：$P'=(TP)_{1:3}=RP+t=\begin{bmatrix}X'\\Y'\\Z'\end{bmatrix}$

$s\hat{u}=KP'$

展开：

$\begin{bmatrix}su\\sv\\s\end{bmatrix}=\begin{bmatrix}f_x&0&c_x\\0&f_y&c_y\\0&0&1\end{bmatrix}\begin{bmatrix}X'\\Y'\\Z'\end{bmatrix}$

即：

$\left\{\begin{matrix}su=f_xX'+c_xZ'\\sv=f_yY'+c_yZ'\\s=Z'\end{matrix}\right.$

整理得：

$\left\{\begin{matrix}u=\frac{f_xX'}{Z'}+c_x \\ v=\frac{f_yY'}{Z'}+c_y\end{matrix}\right.$

### 求解误差项对待估计变量的偏导数
#### 优化位姿
使用扰动模型：
$\frac{\partial e}{\partial \delta \xi}=\mathop{\lim}\limits_{\delta \xi \to 0}\frac{e(\delta \xi \bigoplus \xi)-e(\xi)}{\delta \xi}=\frac{\partial e}{\partial \hat{u}}\frac{\partial \hat{u}}{\partial \delta \xi}=\frac{\partial e}{\partial \hat{u}}\frac{\partial \hat{u}}{\partial P'}\frac{\partial P'}{\partial \delta \xi}$

$\frac{\partial e}{\partial \hat{u}}=\frac{\partial (u-\hat{u})}{\partial \hat{u}}=-I$

$\frac{\partial \hat{u}}{\partial P'}=\begin{bmatrix}\frac{\delta u}{\delta X'}&\frac{\delta u}{\delta Y'}&\frac{\delta u}{\delta Z'}\\\frac{\delta v}{\delta X'}&\frac{\delta v}{\delta Y'}&\frac{\delta v}{\delta Z'}\end{bmatrix}=\begin{bmatrix}\frac{f_x}{Z'}&0&-\frac{f_xX'}{Z'^2}\\0&\frac{f_y}{Z'}&-\frac{f_yY'}{Z'^2}\end{bmatrix}$

$\frac{\partial P'}{\partial \delta \xi}=(\frac{\partial (TP)}{\partial \delta \xi})_{1:3}$

$\frac{\partial (TP)}{\partial \delta \xi}=\mathop{\lim}\limits_{\delta \xi \to 0}\frac{\exp(\delta \xi^\wedge)\exp(\xi^\wedge)P-\exp(\xi^\wedge)P}{\delta \xi}\\=\mathop{\lim}\limits_{\delta \xi \to 0}\frac{(I+\delta \xi^\wedge)\exp(\xi ^ \wedge)P-\exp(\xi^\wedge)P}{\delta \xi}\\=\mathop{\lim}\limits_{\delta \xi \to 0}\frac{\delta \xi^\wedge\exp(\xi ^ \wedge)P}{\delta \xi}\\=\mathop{\lim}\limits_{\delta \xi \to 0}\frac{\begin{bmatrix}\delta\phi^\wedge&\delta\rho\\0^T&0\end{bmatrix}\begin{bmatrix}R&t\\0^T&1\end{bmatrix}P}{\delta \xi}\\=\mathop{\lim}\limits_{\delta \xi \to 0}\frac{\begin{bmatrix}\delta\phi^\wedge&\delta\rho\\0^T&0\end{bmatrix}\begin{bmatrix}RP+t\\1\end{bmatrix}}{\delta \xi}\\=\mathop{\lim}\limits_{\delta \xi \to 0}\frac{\begin{bmatrix}\delta\phi^\wedge(RP+t)+\delta\rho\\0\end{bmatrix}}{\begin{bmatrix}\delta \rho\\\delta \phi\end{bmatrix}}\\=\begin{bmatrix}\frac{\partial (\delta\phi^\wedge(RP+t)+\delta\rho)}{\partial \delta \rho}&\frac{\partial(\delta\phi^\wedge(RP+t)+\delta\rho)}{\partial \delta \phi}\\\frac{\partial 0}{\partial \delta \rho}&\frac{\partial 0}{\partial \delta \phi}\end{bmatrix}\\=\begin{bmatrix}\frac{\partial (\delta\phi^\wedge-(RP+t)+\delta\rho)}{\partial \delta \rho}&\frac{\partial(-(RP+t)^\wedge\delta\phi+\delta\rho)}{\partial \delta \phi}\\\frac{\partial 0}{\partial \delta \rho}&\frac{\partial 0}{\partial \delta \phi}\end{bmatrix}\\=\begin{bmatrix}I&-(RP+t)^\wedge\\0^T&0^T\end{bmatrix}\\\stackrel{\mathrm{def}}{=}(TP)^\odot$

则：

$\frac{\partial P'}{\partial \delta \xi}=(\frac{\partial (TP)}{\partial \delta \xi})_{1:3}=((Tp)^\odot)_{1:3}=\begin{bmatrix}I&-(RP+t)^\wedge\end{bmatrix}=\begin{bmatrix}I&-P'^\wedge\end{bmatrix}=\begin{bmatrix}1&0&0&0&Z'&-Y'\\0&1&0&-Z'&0&X'\\0&0&1&Y'&-X'&0\end{bmatrix}$

$\frac{\partial e}{\partial \delta \xi}=\frac{\partial e}{\partial \hat{u}}\frac{\partial \hat{u}}{\partial P'}\frac{\partial P'}{\partial \delta \xi}\\=-I\begin{bmatrix}\frac{f_x}{Z'}&0&-\frac{f_xX'}{Z'^2}\\0&\frac{f_y}{Z'}&-\frac{f_yY'}{Z'^2}\end{bmatrix}\begin{bmatrix}1&0&0&0&Z'&-Y'\\0&1&0&-Z'&0&X'\\0&0&1&Y'&-X'&0\end{bmatrix}\\=-\begin{bmatrix}\frac{f_x}{Z'}&0&-\frac{f_xX'}{Z'^2}&-\frac{f_xX'Y'}{Z'^2}&f_x+\frac{f_xX'^2}{Z'^2}&-\frac{f_xY'}{Z'}\\0&\frac{f_y}{Z'}&-\frac{f_yY'}{Z'^2}&-f_y-\frac{f_yY'^2}{Z'^2}&\frac{f_yX'Y'}{Z'^2}&\frac{f_yX'}{Z'}\end{bmatrix}$

注1：在李群和李代数的理论中，$\exp(\delta \xi^\wedge)$ 和 $(I+\delta \xi^\wedge)$ 是在某些情况下近似相等的，而不是严格相等。

这个近似来自于泰勒级数的线性近似。泰勒级数是一个无穷级数，它可以用来近似复杂的函数。对于指数函数，它的泰勒级数是：

$$\exp(x) = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \frac{x^4}{4!} + \cdots$$

当 $x$ 很小的时候，高阶项（如 $x^2, x^3, x^4, \cdots$）的贡献会变得非常小，所以我们可以忽略它们，只保留前两项 $1 + x$，这就是线性近似。

所以，当 $\delta \xi^\wedge$ 很小的时候，我们有 $\exp(\delta \xi^\wedge) \approx I+\delta \xi^\wedge$。

但是请注意，这只是一个近似，只有当 $\delta \xi^\wedge$ 很小的时候才成立。当 $\delta \xi^\wedge$ 不小的时候，这个近似可能会有较大的误差。

注2：$a^\wedge b=a\times b=-b \times a=-b^\wedge a$

#### 优化点
$\frac{\partial e}{\partial P}=\frac{\partial e}{\partial \hat{u}}\frac{\partial \hat{u}}{\partial P'}\frac{\partial P'}{\partial P}\\=\frac{\partial e}{\partial \hat{u}}\frac{\partial \hat{u}}{\partial P'}\frac{\partial (RP+t)}{\partial P}=-\begin{bmatrix}\frac{f_x}{Z'}&0&-\frac{f_xX'}{Z'^2}\\0&\frac{f_y}{Z'}&-\frac{f_yY'}{Z'^2}\end{bmatrix}R$
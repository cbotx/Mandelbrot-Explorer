# BigFixed 乘法实现与优化分析

本文档描述当前 `Mandelbrot-Explorer` 中 BigFixed 的实际使用条件、固定点表示、普通乘法、高半乘法、复数平方参考轨道，以及后续可尝试的优化方向。

对应代码：

- `src/engine/bigfixed.h`
- `src/engine/mandel_perturbation.cpp`
- `src/engine/mandel_perturbation.h`
- `src/bench/bench_bigfixed.cpp`
- `src/bench/bench_mulhigh.cpp`

## 1. BigFixed 是否已被实际采用

BigFixed 已经接入生产参考轨道，但只用于 Mandelbrot 或兼容二次公式的深层扰动引擎，不用于通用公式的 MPFR renderer。

当前调度逻辑：

```cpp
if (MANDEL_BIGFIXED 未设置)
    BigFixed = _use_floatexp;
else
    BigFixed = atoi(MANDEL_BIGFIXED) > 0;
```

也就是：

- `MANDEL_BIGFIXED=0`：强制使用 GMP `mpf_t`。
- `MANDEL_BIGFIXED=1`：强制使用 BigFixed。
- 环境变量未设置：FloatExp 深路径自动使用 BigFixed。

FloatExp 通常在以下情况启用：

- zoom 大于约 `1e280`。
- zoom 没有达到 `1e280`，但参考轨道经过非常接近零的位置，普通 double 扰动项可能进入 denormal/underflow 区域。

周期参考轨道仍然使用 `mpf_t`。BigFixed 主要优化完整的串行参考轨道：

```text
z[n+1] = z[n]^2 + c
```

当前真实整帧收益小于独立乘法微基准，因为参考构建只占整帧的一部分。历史 1e876 测量中，参考构建约从 `1.261 s` 降至 `1.169 s`，约为 `1.08x`。

## 2. 固定点表示

定义基数：

\[
B=2^{64}
\]

BigFixed 保存：

```cpp
int sign;
int L;
uint64_t m[L];
```

其中 `m` 是 little-endian limb：

```text
m[0]       最低位 limb
m[L - 1]   最高位 limb
```

对应数学值：

\[
x=\operatorname{sign}(x)\left(\sum_{i=0}^{L-1}m_iB^i\right)B^{-(L-1)}
\]

展开后：

\[
x=\operatorname{sign}(x)\left(
m_{L-1}
\frac{m_{L-2}}{B}
\frac{m_{L-3}}{B^2}
\cdots
\frac{m_0}{B^{L-1}}
\right)
\]

可以把布局理解为：

```text
m[L-1] . m[L-2] m[L-3] ... m[1] m[0]
整数 limb          小数 limbs
```

### 2.1 固定精度

小数精度为：

\[
64(L-1)\text{ bit}
\]

一个 BigFixed ULP 为：

\[
\operatorname{ULP}=B^{-(L-1)}=2^{-64(L-1)}
\]

### 2.2 整数范围

当前只有一个整数 limb，因此正常范围约为：

\[
|x|<2^{64}
\]

代码假设 Mandelbrot 参考轨道的 `z` 和 `c` 都是 `O(1)`，因此加法最高 carry 和乘法最高溢出 limb 可以忽略。

标准 Mandelbrot bailout 为 100，这个假设成立。Custom 二次公式如果允许远大于 `2^64` 的 bailout，则该假设可能失效，后文会单独讨论。

### 2.3 L 的选择

生产代码：

```cpp
int prec = mpf_get_prec(_ref_z_re);
_bfL = (prec + 63) / 64 + 2;
```

含义：

```text
ceil(prec / 64) 个目标精度 limb
+ 1 个整数 limb
+ 1 个额外 guard limb
```

固定点在 `|z|` 非常小时会损失相对精度，因此多保留一个低位 guard limb。

## 3. bf_mul 总调用链

生产入口：

```cpp
bf_mul(result, a, b, scratch);
```

调用关系：

```text
bf_mul
└── bf_mag_mulshift
    ├── L < 16
    │   └── bf_mag_mul_full
    │       └── mpn_mul_n
    └── L >= 16
        └── bf_mag_mulhigh
            └── L 次 mpn_addmul_1
```

`scratch` 由调用者预先分配为 `2L` 个 limb，参考轨道热循环中没有重复分配。

## 4. 乘法的缩放公式

这一节从每个符号的含义开始推导，避免把“真实数值”“无符号 limb 整数”和“固定缩放因子”混在一起。

### 4.1 基本符号

定义：

\[
B=2^{64}
\]

`B` 是 limb 基数。一个 `uint64_t` limb 保存 `0` 到 `B-1`。

定义：

\[
L=\text{每个 BigFixed 的 limb 数量}
\]

定义统一缩放因子：

\[
S=B^{L-1}=2^{64(L-1)}
\]

`S` 决定二进制点的位置。所有具有相同 `L` 的 BigFixed 都除以同一个 `S`。

现在考虑两个输入 BigFixed，分别叫作 `x` 和 `y`。

`x` 的代码字段：

```cpp
x.sign
x.m[0 ... L-1]
```

定义：

\[
s_x=x.\text{sign}\in\{-1,0,+1\}
\]

含义：

- `s_x=+1`：`x` 是正数。
- `s_x=-1`：`x` 是负数。
- `s_x=0`：`x` 是精确零。

定义 `x` 的无符号 limb 整数：

\[
X=\sum_{i=0}^{L-1}x_iB^i
\]

其中：

\[
x_i=x.m[i]
\]

所以 `X` 是把 `x.m[]` 当作一个 little-endian、L-limb、无符号大整数后的结果。

真实数值 `x` 为：

\[
x=s_x\frac{X}{S}
\]

同样，对于输入 `y`：

\[
s_y=y.\text{sign}\in\{-1,0,+1\}
\]

\[
Y=\sum_{i=0}^{L-1}y_iB^i
\]

\[
y=s_y\frac{Y}{S}
\]

代码字段与数学符号的对应关系：

| 代码 | 数学符号 | 含义 |
|---|---:|---|
| `x.sign` | \(s_x\) | `x` 的符号 |
| `y.sign` | \(s_y\) | `y` 的符号 |
| `x.m[i]` | \(x_i\) | `x` 的第 i 个 limb |
| `y.m[i]` | \(y_i\) | `y` 的第 i 个 limb |
| `x.m[]` 组成的大整数 | \(X\) | `x` 的无符号整数 mantissa |
| `y.m[]` 组成的大整数 | \(Y\) | `y` 的无符号整数 mantissa |
| \(B^{L-1}\) | \(S\) | 所有 BigFixed 共用的缩放因子 |

### 4.2 为什么输入要除以 S

以 `L=3` 为例：

\[
S=B^2
\]

无符号整数：

\[
X=x_0+x_1B+x_2B^2
\]

真实数值：

\[
x=s_x\frac{x_0+x_1B+x_2B^2}{B^2}
\]

展开：

\[
x=s_x\left(
x_2+\frac{x_1}{B}+\frac{x_0}{B^2}
\right)
\]

所以：

```text
x.m[2]  是整数 limb
x.m[1]  是第一个小数 limb
x.m[0]  是第二个小数 limb
```

### 4.3 精确乘积

两个真实输入为：

\[
x=s_x\frac{X}{S}
\]

\[
y=s_y\frac{Y}{S}
\]

因此精确乘积：

\[
xy=
\left(s_x\frac{X}{S}\right)
\left(s_y\frac{Y}{S}\right)
\]

\[
xy=s_xs_y\frac{XY}{S^2}
\]

乘积符号为：

\[
s_z=s_xs_y
\]

代码中：

```cpp
r.sign = (x.sign == y.sign) ? +1 : -1;
```

如果任一输入的 sign 为 0，则：

\[
xy=0
\]

函数直接返回精确零。

### 4.4 输出必须继续使用相同的缩放因子

BigFixed 输出 `z` 仍然必须使用相同的格式：

\[
z=s_z\frac{Z}{S}
\]

其中 `Z` 必须是一个最多 L-limb 的无符号整数。

我们希望：

\[
z\approx xy
\]

代入：

\[
s_z\frac{Z}{S}
\approx
s_xs_y\frac{XY}{S^2}
\]

因为：

\[
s_z=s_xs_y
\]

可以去掉共同符号：

\[
\frac{Z}{S}
\approx
\frac{XY}{S^2}
\]

两边乘以 `S`：

\[
Z\approx\frac{XY}{S}
\]

因此输出整数 mantissa 应该是：

\[
Z=\operatorname{round}\left(\frac{XY}{S}\right)
\]

又因为：

\[
S=B^{L-1}
\]

所以：

\[
Z=\operatorname{round}\left(
\frac{XY}{B^{L-1}}
\right)
\]

这就是为什么完整 `2L`-limb 乘积完成后，必须右移：

\[
L-1\text{ 个 limb}
\]

而不是右移 `2(L-1)` 个 limb。

输出在被解释成真实数值时还会再自动除以一次 `S`：

\[
z=s_z\frac{Z}{S}
\]

将：

\[
Z\approx\frac{XY}{S}
\]

代回：

\[
z\approx
s_z\frac{XY/S}{S}
=s_z\frac{XY}{S^2}
=xy
\]

### 4.5 用整数商和余数表示

先计算完整无符号整数积：

\[
P=XY
\]

`P` 最多需要 `2L` 个 limb。

将 `P` 除以缩放因子：

\[
P=QS+T
\]

其中：

\[
Q=\left\lfloor\frac{P}{S}\right\rfloor
\]

\[
T=P\bmod S
\]

范围：

\[
0\le T<S
\]

`Q` 是直接截断后的输出 mantissa，`T` 是被丢弃的低位部分。

round-half-up 的结果为：

\[
Z=
\begin{cases}
Q,&T<\frac{S}{2}\\
Q+1,&T\ge\frac{S}{2}
\end{cases}
\]

### 4.6 为什么只检查一个 discarded bit

因为：

\[
S=B^{L-1}=2^{64(L-1)}
\]

所以：

\[
\frac{S}{2}=2^{64(L-1)-1}
\]

余数 `T` 正好由完整乘积的最低 `L-1` 个 limb 组成：

```text
P[0 ... L-2]
```

其中 `P[L-2]` 的 bit 63 就是：

\[
2^{64(L-1)-1}
\]

因此：

```cpp
P[L - 2] >> 63
```

直接回答：

```text
T 是否 >= S/2
```

当前舍入规则是 ties-up，所以：

- 最高 discarded bit 为 0：一定小于半个 ULP。
- 最高 discarded bit 为 1：一定大于或等于半个 ULP。

不需要 sticky bit 区分“刚好一半”和“大于一半”，因为两者都会向上。

### 4.7 最终误差

舍入保证：

\[
\left|
Z-\frac{XY}{S}
\right|
\le\frac12
\]

输出真实值：

\[
z=s_z\frac{Z}{S}
\]

精确乘积：

\[
xy=s_z\frac{XY/S}{S}
\]

所以：

\[
|z-xy|
=
\frac{1}{S}
\left|
Z-\frac{XY}{S}
\right|
\le
\frac{1}{2S}
\]

而一个 BigFixed ULP 是：

\[
\operatorname{ULP}=\frac1S
\]

因此完整正确舍入乘法的绝对误差最多为：

\[
0.5\operatorname{ULP}
\]

### 4.8 十进制类比示例

实际代码使用：

\[
B=2^{64}
\]

为了容易观察，这里临时使用十进制类比：

```text
B = 100
L = 3
S = B^(L-1) = 10000
```

表示：

```text
x = +1.25
y = -0.50
```

选择整数 mantissa：

\[
X=12500
\]

\[
Y=5000
\]

符号：

\[
s_x=+1
\]

\[
s_y=-1
\]

检查输入：

\[
x=s_x\frac{X}{S}
=+\frac{12500}{10000}
=1.25
\]

\[
y=s_y\frac{Y}{S}
=-\frac{5000}{10000}
=-0.5
\]

整数乘积：

\[
P=XY=12500\cdot5000=62500000
\]

重新缩放：

\[
Z=\operatorname{round}\left(\frac{P}{S}\right)
=\operatorname{round}\left(\frac{62500000}{10000}\right)
=6250
\]

输出符号：

\[
s_z=s_xs_y=-1
\]

输出真实值：

\[
z=s_z\frac{Z}{S}
=-\frac{6250}{10000}
=-0.625
\]

与精确结果一致：

\[
1.25\cdot(-0.5)=-0.625
\]

### 4.9 与代码步骤一一对应

数学：

\[
P=XY
\]

代码：

```cpp
mpn_mul_n(tmp, x.m, y.m, L);
```

数学：

\[
Q=\left\lfloor\frac{P}{B^{L-1}}\right\rfloor
\]

代码：

```cpp
const int drop = L - 1;
for (int i = 0; i < L; ++i)
    result.m[i] = tmp[drop + i];
```

数学：

\[
Z=
\begin{cases}
Q,&T<S/2\\
Q+1,&T\ge S/2
\end{cases}
\]

代码：

```cpp
if (tmp[drop - 1] >> 63)
    add_one_with_carry(result.m);
```

数学：

\[
s_z=s_xs_y
\]

代码：

```cpp
result.sign =
    (x.sign == y.sign) ? +1 : -1;
```

总结后的核心步骤：

```text
输入：
x = s_x * X / S
y = s_y * Y / S

整数乘积：
P = X * Y

重新缩放：
Z = round(P / S)

输出：
z = (s_x * s_y) * Z / S
```

## 5. 小 L：完整乘积路径

当：

```cpp
L < 16
```

使用 GMP 的底层整数乘法：

```cpp
mpn_mul_n(tmp, a, b, L);
```

输出是完整的 `2L` limb：

```text
tmp[0], tmp[1], ..., tmp[2L-1]
```

### 5.1 计算完整整数积

\[
P=A\cdot C
\]

其中：

\[
P=\sum_{k=0}^{2L-1}p_kB^k
\]

### 5.2 固定点右移

需要右移：

\[
F=L-1
\]

个 limb。

代码：

```cpp
const int drop = L - 1;
for (int i = 0; i < L; ++i)
    r[i] = tmp[drop + i];
```

布局：

```text
丢弃：
tmp[0 ... L-2]

保留：
tmp[L-1 ... 2L-2]
```

`tmp[2L-1]` 是超出当前一个整数 limb 范围的最高溢出部分。当前实现依赖 `|z|=O(1)` 的假设，不保存它。

### 5.3 舍入

保留区最低 limb 是：

```text
tmp[L-1]
```

被丢弃部分的最高 limb 是：

```text
tmp[L-2]
```

检查它的最高 bit：

```cpp
if (tmp[drop - 1] >> 63)
    result += 1;
```

如果该 bit 为 1，被丢弃值至少达到半个 ULP，因此保留结果加一。

进位传播：

```cpp
unsigned char carry = _addcarry_u64(0, r[0], 1, &r[0]);
for (int i = 1; i < L && carry; ++i)
    carry = _addcarry_u64(carry, r[i], 0, &r[i]);
```

当前舍入规则是：

```text
round to nearest, ties upward
```

不是 ties-to-even。

对于这种舍入规则，只检查最高 discarded bit 即可：

- bit 为 0：严格小于 0.5 ULP。
- bit 为 1：大于或等于 0.5 ULP。

不需要 sticky bit 来区分 tie 和大于 tie，因为两者都向上。

## 6. 大 L：High-half 短乘法

当：

```cpp
L >= 16
```

调用：

```cpp
bf_mag_mulhigh(tmp, a, b, L, 3);
```

定义：

```cpp
F = L - 1;
GUARD = 3;
lo = F - GUARD = L - 4;
```

完整乘积需要的列为：

```text
0 ... 2L-1
```

短乘法只显式累加：

```text
lo ... 2L-1
```

也就是：

```text
L-4 ... 2L-1
```

其中：

```text
L-4, L-3, L-2   三个 guard limb
L-1 ... 2L-2    最终 L 个结果 limb
```

### 6.1 hbuf 的列映射

临时数组使用局部下标：

```text
hbuf[t] = 完整乘积的 column (lo + t)
```

因此：

```text
hbuf[0]       -> column L-4
hbuf[1]       -> column L-3
hbuf[2]       -> column L-2
hbuf[3]       -> column L-1
...
```

### 6.2 遍历 a 的每个 limb

```cpp
for (int i = 0; i < L; ++i)
```

第 `i` 行包含：

\[
a_i b_0,\ a_i b_1,\ldots,a_i b_{L-1}
\]

第 `i,j` 个部分积属于完整乘积列：

\[
k=i+j
\]

### 6.3 找到第一个需要计算的 j

只保留：

\[
i+j\ge lo
\]

所以：

```cpp
int jstart = lo - i;
if (jstart < 0)
    jstart = 0;
```

低于 `jstart` 的部分积全部跳过。

当前行需要计算的 limb 数：

```cpp
int n = L - jstart;
```

### 6.4 映射到 hbuf

第一个部分积的完整列为：

\[
i+jstart
\]

对应 hbuf：

```cpp
int t = i + jstart - lo;
```

### 6.5 调用 mpn_addmul_1

```cpp
mp_limb_t carry = mpn_addmul_1(
    hbuf + t,
    b + jstart,
    n,
    a[i]);
```

等价于：

```text
hbuf[t ... t+n-1] += b[jstart ... L-1] * a[i]
```

`mpn_addmul_1` 是 GMP 已优化的单 limb 乘整向量并累加 primitive，Windows x64 版本会使用 GMP 的汇编实现。

### 6.6 传播行末 carry

`mpn_addmul_1` 返回最高 carry：

```cpp
int u = t + n;
while (carry) {
    hbuf[u] += carry;
    carry = carry_out;
    ++u;
}
```

### 6.7 提取 L 个结果 limb

最终保留区从完整列 `L-1` 开始。

因为：

```text
lo = L-4
GUARD = 3
```

所以：

```cpp
for (int i = 0; i < L; ++i)
    r[i] = hbuf[i + GUARD];
```

即：

```text
r[0] = hbuf[3] = column L-1
r[1] = hbuf[4] = column L
...
```

### 6.8 使用 guard limb 舍入

结果最低 limb 下面的第一个完整 limb 是：

```text
hbuf[GUARD-1] = hbuf[2] = column L-2
```

代码：

```cpp
if (hbuf[GUARD - 1] >> 63)
    result += 1;
```

然后通过 `_addcarry_u64` 向高位传播。

## 7. High-half 的严格误差问题

代码注释目前把结果高 L limb 描述为 `EXACT`。随机测试也没有观察到与完整 `mpn_mul_n` 的高位差异，但从严格数学角度看，仅保留有限 guard limb 并不能直接证明结果永远 exact。

原因是被省略的低列：

```text
0 ... lo-1
```

仍然可能产生向 `lo` 列的 carry。这个 carry 继续穿过三个 guard limb 后，理论上可能影响最终最低保留 limb或舍入判定。

设被省略的全部非负部分积之和为 `D`。粗略上界：

\[
D < L B^{lo+1}
\]

目标结果一个 ULP 对应：

\[
B^F
\]

其中：

\[
F=L-1,\qquad lo=F-G
\]

因此：

\[
\frac{D}{B^F}<L B^{1-G}
\]

当前：

\[
G=3
\]

所以：

\[
\frac{D}{B^F}<\frac{L}{B^2}
=\frac{L}{2^{128}}\text{ ULP}
\]

这个误差极小，因此随机测试几乎不可能碰到差异。但只要 truncated product 恰好距离舍入边界小于这个上界，低位 carry 仍可能改变最终一 ULP。

因此更准确的描述是：

```text
三个 guard limb 让错误概率和误差上界极小，但尚未形成严格的 correctly-rounded 证明。
```

## 8. bf_mul 的符号与零

入口：

```cpp
inline void bf_mul(
    BigFixed& r,
    const BigFixed& a,
    const BigFixed& b,
    uint64_t* tmp);
```

步骤：

```cpp
if (a.sign == 0 || b.sign == 0) {
    r = 0;
    return;
}

bf_mag_mulshift(r.m, a.m, b.m, L, tmp);

r.sign = (a.sign == b.sign) ? +1 : -1;
```

如果 magnitude 在固定精度下被舍入成全零，当前代码仍可能保留非零 sign。所有转换函数都会扫描 magnitude，并把全零 magnitude 当成数值零，因此当前行为不会改变数值。

更严格的表示不变量可以在乘法后检查最低到最高 limb 是否全零，再将 `sign=0`，但这会增加每次热乘法的 O(L) 扫描，所以当前没有这样做。

## 9. 生产复数平方

参考轨道：

\[
z_{n+1}=z_n^2+c
\]

设：

\[
z=a+bi
\]

则：

\[
z^2=(a^2-b^2)+2ab\,i
\]

生产代码没有直接计算三次乘法：

```text
a*a
b*b
a*b
```

而是使用两乘法 Karatsuba 形式：

\[
a^2-b^2=(a+b)(a-b)
\]

完整步骤：

```cpp
s  = a + b;
d  = a - b;

ab = round(a * b);
re = round(s * d);

im = ab + ab;

nextRe = re + cRe;
nextIm = im + cIm;
```

对应代码：

```cpp
bf_add(_bt1, a, b);
bf_sub(_bt2, a, b);
bf_mul(_bab, a, b, scratch);
bf_mul(_bre, _bt1, _bt2, scratch);
bf_add(_bim, _bab, _bab);
bf_add(nextRe, _bre, cRe);
bf_add(nextIm, _bim, cIm);
```

每轮成本：

```text
2 次 BigFixed 乘法
5 次 BigFixed 加减
2 次 BigFixed -> FloatExp shadow 提取
```

所有 BigFixed 对象和乘法 scratch 都在参考构建开始前分配，热循环内没有 vector 扩容。

## 10. 舍入顺序

Karatsuba 实部：

\[
(a+b)(a-b)=a^2-b^2
\]

`a+b` 和 `a-b` 是固定点精确加减，然后乘法只舍入一次。

虚部：

\[
2ab
\]

当前先计算并舍入：

\[
\widehat{ab}=\operatorname{round}(ab)
\]

再执行精确固定点加法：

\[
\widehat{\Im}=2\widehat{ab}
\]

这不一定等于：

\[
\operatorname{round}(2ab)
\]

换成三平方公式时，舍入位置也会改变，因此即使最终误差只有 1–2 ULP，混沌轨道后续也可能产生不同 escape iteration。

## 11. BigFixed 与 mpf 的差异

BigFixed：

- 固定绝对精度。
- 乘法 round-half-up。
- 加减在固定网格上精确。

GMP `mpf_t`：

- 浮点格式，具有动态 exponent。
- 当前参考路径的运算舍入行为更接近截断。
- 接近零时仍保持相对精度。

因此 BigFixed 参考轨道不保证和 `mpf_t` byte-identical。验收必须以更高精度 oracle 的最终分类和 escape iteration 为准，而不能只比较 BigFixed 与 `mpf_t` 的内部轨道字节。

## 12. 当前微基准

测试机器：

```text
AMD EPYC 7763
8 cores / 16 threads
MSVC 19.51
GMP static x64
```

### 12.1 单次 raw multiply

| L | 精度 | BigFixed | GMP mpf | BigFixed / mpf |
|---:|---:|---:|---:|---:|
| 3 | 192 bit | 16.5 ns | 9.6 ns | 0.58x |
| 9 | 576 bit | 74.9 ns | 88.6 ns | 1.18x |
| 24 | 1536 bit | 448.7 ns | 481.7 ns | 1.07x |
| 48 | 3072 bit | 1346.9 ns | 1452.0 ns | 1.08x |

低精度 raw multiply 中 BigFixed 可能比 `mpf_mul` 慢。生产自动调度只在深 FloatExp 路径开启，不在浅层启用。

### 12.2 生产一致的两乘法复数轨道步

| L | 精度 | mpf 两乘法 | BigFixed 两乘法 | 加速 |
|---:|---:|---:|---:|---:|
| 3 | 192 bit | 79.1 ns | 42.3 ns | 1.87x |
| 9 | 576 bit | 142.2 ns | 111.9 ns | 1.27x |
| 24 | 1536 bit | 561.5 ns | 400.3 ns | 1.40x |
| 48 | 3072 bit | 1558.3 ns | 1182.8 ns | 1.32x |

即使 raw multiply 加速有限，固定点加减、无 exponent bookkeeping 和更轻量的 shadow 提取仍会改善完整参考轨道步。

### 12.3 High-half 与完整 mpn_mul_n

| L | mpn_mul_n | addmul high-half | high-half 相对完整 |
|---:|---:|---:|---:|
| 16 | 182.7 ns | 175.5 ns | 1.04x |
| 20 | 287.4 ns | 246.9 ns | 1.16x |
| 24 | 423.3 ns | 339.2 ns | 1.25x |
| 32 | 602.1 ns | 591.9 ns | 1.02x |
| 48 | 1376.9 ns | 1107.3 ns | 1.24x |
| 64 | 1896.4 ns | 1883.2 ns | 1.01x |
| 96 | 3569.3 ns | 4064.0 ns | 0.88x |
| 128 | 6151.6 ns | 7150.3 ns | 0.86x |

当前 `L>=16` 始终使用 high-half，但在 L≈96 以后，GMP 完整乘法的 Karatsuba/Toom 已经击败 O(L²) 的逐行 `mpn_addmul_1`。

## 13. 可优化方向

## 13.1 为 high-half 增加严格的舍入认证

目标：

```text
常见情况继续使用短乘法
只有舍入边界不确定时回退完整 mpn_mul_n
```

可行方案：

1. 短乘法计算三个 guard limb和保留高位。
2. 为所有省略部分积计算严格非负上界 `E`。
3. 当前 truncated product 为下界 `Q`。
4. 真正乘积属于：

\[
P\in[Q,Q+E]
\]

5. 检查整个区间是否会跨越：
   - 保留整数边界。
   - 0.5 ULP 舍入边界。
6. 如果 `[Q,Q+E]` 的所有值舍入结果相同，则直接接受。
7. 如果可能产生两个不同结果，则调用完整 `mpn_mul_n`。

当前 `G=3` 时：

\[
E<\frac{L}{2^{128}}\text{ ULP}
\]

因此 fallback 概率应当极低，正常性能几乎不受影响。

这是最优先的正确性改进。

## 13.2 增加 high-half 的上界 crossover

当前：

```cpp
if (L >= 16)
    high_half();
else
    mpn_mul_n();
```

根据当前机器的测量，可以先尝试：

```cpp
if (L >= 16 && L < 80)
    high_half();
else
    mpn_mul_n();
```

更保守的静态候选：

```text
16 <= L < 72   high-half
其他           full mpn_mul_n
```

应在不同机器上重新测量。也可以在程序启动时做一次极短的 kernel calibration，但这会增加调度复杂度和结果复现成本。

## 13.3 实现 high-half square

当前 `bf_sqr` 使用完整：

```cpp
mpn_sqr
```

`bench_mulhigh.cpp` 已有实验性 `sqrhigh_addmul`：

1. 只计算一次上三角非对角项。
2. 将非对角和乘二。
3. 加入对角项 `a[i]^2`。
4. 只保留高半区与 guard limbs。

复数平方可以改成：

\[
s_1=a^2
\]

\[
s_2=b^2
\]

\[
s_3=(a+b)^2
\]

\[
\Re=s_1-s_2
\]

\[
\Im=s_3-s_1-s_2
\]

纯乘法测量：

| L | 两次 high-half mul | 三次 high-half sqr | 三平方相对两乘法 |
|---:|---:|---:|---:|
| 16 | 363.8 ns | 435.8 ns | 0.83x |
| 24 | 692.3 ns | 785.5 ns | 0.88x |
| 32 | 1129.7 ns | 1192.1 ns | 0.95x |
| 48 | 2386.0 ns | 2665.3 ns | 0.90x |
| 64 | 3937.3 ns | 3498.9 ns | 1.13x |

因此只值得从 L≈64 以上测试。

注意三平方和两乘法的舍入顺序不同，不能只做单步 ULP 测试，还必须做完整轨道和最终图像 oracle 验证。

## 13.4 融合 a+b 与 a-b

当前：

```cpp
bf_add(s, a, b);
bf_sub(d, a, b);
```

根据 sign-magnitude 情况，两次调用可能分别执行 magnitude compare 和完整 limb 扫描。

可以实现：

```cpp
bf_sumdiff(s, d, a, b);
```

在一次准备阶段中：

1. 比较一次 `|a|` 与 `|b|`。
2. 根据符号决定 sum/difference 的 magnitude 运算方向。
3. 在一轮或两轮紧凑 limb loop 中同时生成两个结果。

可能减少：

- 一次 magnitude compare。
- 一部分重复读取。
- 函数和分支开销。

预期收益为几个百分点，低于乘法 kernel 优化。

## 13.5 融合两次 high-half multiply

当前每个复数平方分别调用：

```text
mul(a, b)
mul(a+b, a-b)
```

每次调用都会：

- 清零 hbuf。
- 遍历 L 行。
- 调用 L 次 `mpn_addmul_1`。

可以考虑双输出 kernel：

```cpp
bf_mag_mulhigh_pair(
    outAB, a, b,
    outSD, sum, difference);
```

潜在收益：

- 合并外层行循环。
- 减少循环边界和地址计算。
- 改善 a/b/sum/difference 的 cache locality。

但底层仍需要约 `2L` 次 `mpn_addmul_1`，收益可能有限，需要微基准确认。

## 13.6 超大 L 使用截断 Karatsuba/middle product

当前 high-half 是 O(L²) 的 schoolbook 行累加。完整 `mpn_mul_n` 在大 L 使用 Karatsuba、Toom 或更高级算法，因此 L≈96 后重新领先。

理想实现是：

```text
truncated Karatsuba high product
或 middle product
```

只递归计算影响目标高半区的子乘积。

GMP 公共 `gmp.h` 没有暴露 `mpn_mulhigh_n` 或 `mpn_mulmid`，因此需要：

- 自己实现截断 Karatsuba；或
- 使用 GMP 非公开内部接口，代价是版本耦合和可移植性下降。

这是潜在收益最大、实现复杂度也最高的方向。

## 13.7 Custom 大 bailout 的范围保护

当前一个整数 limb 的表示只能安全覆盖：

\[
|x|<2^{64}
\]

但 Custom 二次公式允许的 bailout 可能远大于该范围。

风险：

- `bf_mag_add` 顶部 carry 被忽略。
- `bf_mag_mulshift` 不保存 `tmp[2L-1]`。
- 数值可能静默 wrap/truncate。

建议至少加入调度门：

```cpp
if (customBailout >= safeBigFixedBound)
    disable BigFixed;
```

更完整的方案是支持多个整数 limb，并将固定小数点位置从 `L-1` 改成独立参数：

```cpp
fractionLimbs
integerLimbs
```

但这会让当前简单的统一右移逻辑复杂化。

## 13.8 不建议继续的方向

### 手写 BMI2/ADX schoolbook

当前微基准：

```text
ADX high-half 只有 mpf 的约 0.43x–0.63x
```

明显慢于 GMP 的 `mpn_addmul_1`，当前实现不值得采用。

### OpenMP 并行单次参考乘法

参考轨道迭代存在严格串行依赖：

```text
z[n+1] 依赖 z[n]
```

只能并行 limb 乘法。当前精度下一次乘法非常短，OpenMP barrier 成本约为乘法本身的 17–23 倍，因此不可行。

### 所有 L 都使用三次完整 mpn_sqr

当前 engine-faithful 测量中，L≤48 的三平方均慢于两乘法 Karatsuba，不能无条件替换。

## 14. 推荐实验顺序

### 第一阶段：正确性与低风险调度

1. 为 high-half 推导并实现严格 omitted-product 上界。
2. 舍入不确定时回退完整 `mpn_mul_n`。
3. 增加 L 上界 crossover，例如 `L<80` 才使用 high-half。
4. 增加 Custom bailout 的 BigFixed 安全范围 gate。

### 第二阶段：中等风险算术优化

1. 实现 `bf_sumdiff`。
2. 测试双输出 high-half multiply。
3. 实现 high-half square。
4. 只在 L≥64 测试三平方复数公式。

### 第三阶段：大精度算法

1. 实现 truncated Karatsuba。
2. 与 `mpn_mul_n` 在 L=64、96、128、192、256 对比。
3. 为不同 L 建立静态 kernel crossover 表。

## 15. 每项优化的验收条件

单步算术：

- 与完整 `mpn_mul_n` 的 correctly-rounded 结果一致，或误差有严格上界。
- 正负号、精确零、接近舍入边界和最高 carry 有专门测试。
- 随机测试之外增加构造出的 carry-chain 极端输入。

完整参考轨道：

- shallow、deep51、ticktock、flake。
- exterior1000、parity1000。
- deep876、minibrot875、extref875。
- 普通 Smooth、EDE、Normal、Feather、Orbit Trap。

最终图像：

- class mismatch 必须为 0。
- 不允许增加 EMPTYPIXEL。
- escape iteration 差异必须在现有 oracle gate 内。
- BLA skip 率不能显著下降。
- 参考构建加速必须转化成整帧收益，不能只看 raw multiply。

## 16. 当前结论

BigFixed 已经被实际采用，而且在 FloatExp 深参考构建中有真实收益。当前两乘法 Karatsuba 复数平方在约 192–3072 bit 范围内优于 `mpf_t`，是合理的生产默认。

最值得优先处理的并不是再写一版普通 schoolbook 乘法，而是：

1. 给 high-half 增加严格舍入认证。
2. 在 L≈80–96 以后切回 `mpn_mul_n`。
3. 保护一个整数 limb 无法覆盖的超大 Custom bailout。
4. 在 L≥64 测试 high-half square 与三平方复数公式。

这些方向分别覆盖正确性、超大精度性能和生产适用范围，风险和收益都比继续微调普通 `_umul128` 循环更明确。

## 17. MSVC 生成汇编检查

本节检查正式构建使用的 MSVC `/O2 /arch:AVX2 /MT` 代码生成，而不是只分析 C++ 源码。

生成汇编清单的命令：

```bat
cmd /v:on /c "call scripts\msvcenv.bat >nul && cl /nologo /O2 /EHsc /openmp /std:c++17 /arch:AVX2 /MT /I \"!VCPKG!\include\" /I src\engine /I src\gui /c src\engine\mandel_perturbation.cpp /Fobuild\asm_inspect\mandel_perturbation.obj /FAs /Fabuild\asm_inspect\mandel_perturbation.asm"
```

微基准：

```bat
cmd /v:on /c "call scripts\msvcenv.bat >nul && cl /nologo /O2 /EHsc /std:c++17 /arch:AVX2 /MT /I \"!VCPKG!\include\" /I src\engine /c src\bench\bench_bigfixed.cpp /Fobuild\asm_inspect\bench_bigfixed.obj /FAs /Fabuild\asm_inspect\bench_bigfixed.asm"
```

### 17.1 完整乘法包装已经是零开销 tail call

`bf_mag_mul_full` 编译为：

```asm
bf_mag_mul_full:
    jmp __gmpn_mul_n
```

说明：

- 没有额外 prologue/epilogue。
- 没有参数搬运。
- 直接 tail-call GMP。
- 完整乘法路径的包装函数不存在可优化开销。

真正的 limb 乘法位于静态 `gmp.lib` 中，由 GMP 汇编实现。

### 17.2 High-half 外层循环

关键结构：

```asm
call    memset

loop_i:
    mov     r9, [a + i*8]

    mov     eax, lo
    sub     eax, i
    xor     edx, edx
    cmovns  edx, eax

    ; edx = max(0, lo-i)
    ; r8d = L-jstart
    ; rcx = hbuf+t
    ; rdx = b+jstart
    ; r9  = a[i]

    call    __gmpn_addmul_1

    test    rax, rax
    je      no_carry

carry_loop:
    add     [hbuf + u*8], rax
    setb    al
    inc     u
    movzx   eax, al
    test    rax, rax
    jne     carry_loop

no_carry:
    inc     i
    cmp     i, L
    jl      loop_i
```

优点：

- `max(0, lo-i)` 编译成 `cmovns`，没有难预测分支。
- `i`、`L`、`lo`、`nout`、数组基址都保存在寄存器中。
- 主乘加由 GMP 汇编 `mpn_addmul_1` 完成。
- carry 通常为零或只传播一个 limb，循环短。
- 没有热循环内存分配。

不足：

- 每一行都调用一次 `mpn_addmul_1`，每次 high-half 乘法有 L 次函数调用。
- 两次复数乘法合计约 `2L` 次 `mpn_addmul_1` 调用。
- 每次乘法先调用一次 `memset` 清空临时缓冲区。
- `mpn_addmul_1` 调用阻止编译器跨行调度和软件流水。

### 17.3 结果复制

源码：

```cpp
for (int i = 0; i < L; ++i)
    r[i] = tmp[i + GUARD];
```

MSVC 生成两条路径：

```asm
; 检查 r 与 tmp 区间是否可能重叠
cmp ...
ja  scalar_copy
cmp ...
jae scalar_copy

call memcpy
```

编译器不知道 `r` 与 `tmp` 必然不重叠，因此生成：

- 两次地址范围比较。
- `memcpy` 快路径。
- 标量复制 fallback。

生产调用中两块内存始终独立，可以直接写：

```cpp
std::memcpy(r, tmp + GUARD, sizeof(uint64_t) * L);
```

或者给参数添加有效的 `__restrict` 契约，从而删除运行时 alias 检查。该优化收益很小，但风险低。

### 17.4 舍入进位

最高 discarded bit 检查编译为：

```asm
cmp     qword ptr [tmp + 16], 0
jge     no_round
```

这是利用有符号判断读取 bit 63：

- 非负：最高 bit 为 0。
- 负：最高 bit 为 1。

进位：

```asm
add     qword ptr [r], 1
setb    carry

round_loop:
    test    carry, carry
    je      done
    adc     limb, 0
    setb    carry
```

这是合理的串行 carry chain，没有多余大整数操作。

### 17.5 bf_mul 包装层

`bf_mul` 没有内联进生产 `calCoefficient`。

每次调用先执行：

```asm
; r.L = L
mov [r], L

; 检查 r.m.size() == L
mov rax, [r.end]
sub rax, [r.begin]
sar rax, 3
cmp rax, L
je  already_sized

call vector::resize
```

生产对象已经在参考构建开始前 `setL(L)`，因此：

```text
vector::resize 分支永远不走
```

但每次乘法仍会执行：

- 读取 vector begin/end。
- 指针相减。
- 右移除以 8。
- 与 L 比较。
- 条件跳转。

然后检查两个 sign：

```asm
cmp [a.sign], 0
je  zero
cmp [b.sign], 0
je  zero
```

最后：

```asm
call bf_mag_mulshift
```

对于 L=24 或 L=48，乘法本身较重，包装开销占比不大；对于 L≈9–16，包装和函数调用占比会明显上升。

可以增加仅供预分配热循环使用的入口：

```cpp
bf_mul_ready(
    uint64_t* result,
    int& resultSign,
    const uint64_t* a,
    int aSign,
    const uint64_t* b,
    int bSign,
    int L,
    uint64_t* scratch);
```

它可以删除：

- vector size 检查。
- `r.L = L` 写入。
- 潜在 resize 分支。
- BigFixed 对象字段重复解引用。

### 17.6 bf_addsub 汇编

同号加法核心编译为：

```asm
loop:
    mov     limb, [a + i*8]
    add     carry_byte, -1
    adc     limb, [b + i*8]
    mov     [r + i*8], limb
    setb    carry_byte
    inc     i
    cmp     i, L
    jl      loop
```

异号减法使用：

```asm
sbb
setb
```

这是正常的串行 carry/borrow chain，单条 limb 流程已经较紧凑。

主要浪费不在单条 `adc/sbb`，而在：

- `bf_addsub` 未内联。
- 每次调用都有 6 个非易失寄存器 push/pop。
- 每次调用检查 vector size。
- 异号时先从最高 limb 做 magnitude compare，再从最低 limb 做 subtraction，需要两次遍历。

生产每轮参考轨道会调用 `bf_addsub` 五次：

```text
a+b
a-b
2ab
real+c.real
imag+c.imag
```

因此 `bf_sumdiff(a,b)` 和预分配 raw-limb API 比继续微调 `adc` 更值得尝试。

### 17.7 calCoefficient 中的调用数量

BigFixed 分支没有被整体内联。

每次参考轨道迭代调用：

```text
2 × bf_addsub     生成 a+b、a-b
2 × bf_mul        生成 ab、(a+b)(a-b)
3 × bf_addsub     生成 2ab、加 c.real、加 c.imag
2 × bf_to_fe      生成 real/imag FloatExp shadow
```

总计：

```text
9 个 BigFixed helper 调用
```

此外：

```text
2 × bf_mul -> 2 × bf_mag_mulshift
2 × bf_mag_mulshift -> 2L × mpn_addmul_1
```

所以生产热步的调用层级为：

```text
calCoefficient
├── bf_addsub
├── bf_addsub
├── bf_mul
│   └── bf_mag_mulshift
│       └── L × mpn_addmul_1
├── bf_mul
│   └── bf_mag_mulshift
│       └── L × mpn_addmul_1
├── bf_addsub
├── bf_addsub
├── bf_addsub
├── bf_to_fe
└── bf_to_fe
```

### 17.8 bf_to_fe 汇编

每次转换先从最高 limb 向下扫描：

```asm
scan:
    mov value, [m + top*8]
    test value, value
    jne found
    dec top
    jns scan
```

然后取最高两个非零附近 limb。

`uint64_t -> double` 对最高 limb 和次高 limb 各有一套 unsigned 转换分支：

```asm
test    limb, limb
js      unsigned_high_bit_path
vcvtsi2sd ...
```

最后调用：

```asm
call frexp
```

生产每轮调用两次 `bf_to_fe`，因此也调用两次 `frexp`。

这部分不是大整数乘法，但在中低 L 时可能占明显比例。

可直接从最高 limb 的 leading-zero count 构造 FloatExp：

1. `_BitScanReverse64` 或 `lzcnt` 得到最高置位 bit。
2. 从 `hi/lo` 拼出最高 53 位。
3. 直接构造 `[0.5,1)` double mantissa。
4. exponent 由 `top`、最高 bit 位置和固定二进制点直接计算。

这样可以删除：

- 两次 `frexp` 调用。
- 大部分 unsigned integer-to-double 分支。
- `hi*2^64+lo` 的大范围浮点归一化。

这是汇编检查后最值得单独微基准的方向之一。

### 17.9 calCoefficient 本身

`calCoefficient` prologue：

```asm
push rbp
push rsi
push rdi
push r12
push r13
push r14
push r15
sub  rsp, 576
```

该函数同时包含：

- BigFixed reference。
- GMP reference。
- FloatExp shadow。
- double shadow。
- derivative/EDE。
- escape continuation。

因此栈帧较大，寄存器压力较高。

每次迭代都在函数内部检查：

```asm
cmp byte ptr [this + use_bigfixed], 0
je  mpf_path
```

可以考虑在参考轨道外层分派一次：

```cpp
if (_use_bigfixed)
    buildReferenceBigFixed();
else
    buildReferenceMpf();
```

然后让 BigFixed 热循环调用更小的专用 step 函数。

潜在收益：

- 删除每迭代 backend 分支。
- 减少 `calCoefficient` 栈帧。
- 让编译器看到 BigFixed 对象均已预分配。
- 更容易内联 `bf_addsub` 包装并消除 size 检查。

但这属于较大重构，应该排在 `bf_to_fe` 和预分配 raw API 的微基准之后。

## 18. 汇编效率结论

### 已经较好的部分

- 完整乘法包装是零开销 tail call。
- high-half 的 `jstart=max(0,lo-i)` 是 branchless `cmov`。
- limb 加减正确使用 `adc/sbb`。
- 舍入 carry chain 紧凑。
- 主 limb 乘加调用 GMP 汇编 `mpn_addmul_1`。
- 没有热循环动态分配。

### 仍有明显开销的部分

- `calCoefficient` 每轮有 9 个 BigFixed helper 调用。
- `bf_mul` 和 `bf_addsub` 每轮重复检查 vector size。
- `bf_addsub` 有较重的非易失寄存器保存/恢复。
- 两次乘法分别清零 scratch。
- 结果复制存在运行时 alias 检查。
- 每轮两次 `bf_to_fe`，每次都扫描 limb 并调用 `frexp`。
- 大 L 的 high-half 仍是 O(L²)，无法利用 `mpn_mul_n` 的 Karatsuba/Toom。

### 基于汇编的新优化优先级

1. 为 `bf_to_fe` 实现基于 `lzcnt/BitScanReverse64` 的无 `frexp` 路径。
2. 增加预分配专用 `bf_mul_ready/bf_addsub_ready`，删除 vector size 检查。
3. 用 `bf_sumdiff` 融合 `a+b` 与 `a-b`。
4. 显式 `memcpy` 或 restrict 契约，删除结果复制 alias 检查。
5. 将 BigFixed reference loop 从通用 `calCoefficient` 中拆出，减少调用和栈帧。
6. 再处理 high-half 严格舍入、L 上界 crossover 和 high-half square。

汇编中没有发现编译器把核心 limb 循环生成为明显低效的代码。当前主要问题不是某一条 `mul/adc` 指令选错，而是 C++ 对象包装、函数边界、重复检查和 shadow 转换叠加出的调度开销。

---
layout: post
title: "표본분산이 chisquare(n-1)을 따르는 이유"
categories: [stat]
tags: [표본분산, 카이제곱, 자유도]
description: 
toc:
  sidebar: left
  beginning: true
---

정규모집단에서 $n$개를 추출한 표본에서 표본분산 $ S^2 = \sum_{i=1}^n ({X_i - \bar{X})} $ 일 때 $ \frac{(n-1)S^2}{\sigma^2} \sim \chi^2 (n-1) $ 라는 것을 배웠을 것이다. 하지만 이것이 왜 진짜로 성립하는지는 입문 수준에서는 잘 다루지 않는다. (보통 저 식을 전개를 해서, $ \chi^2(n) - \chi^2(1) = \chi^2(n-1) $ 이 된다고 설명한다. 이것도 $\chi^2(n)$과 $\chi^2(1)$이 독립일 때 성립하는 것인데, 입문 수준에서는 독립인 이유를 잘 설명하지 않는다.)

아무튼 오늘은 선형대수학을 이용해서 저 사실을 엄밀하게 증명해보려고 한다. 이 글은 독자가 선형대수학에 대한 배경지식이 있을 것이라 가정한다.

<br>
<br>

### 1. 표본의 벡터형 표현

먼저 벡터와 행렬을 이용한 표본과 분포의 표현을 알고 가자. 오늘은 표본분산을 다루는데 집중할 것이므로 이 상황에 맞춰 설명하겠다. 정규모집단 $ N(\mu, \sigma^2) $ 에서 $n$개의 독립인 표본을 뽑았을 때, 우리는 $ X_i \overset{\text{i.i.d.}}{\sim} N(\mu, \sigma^2) $ 라고 썼다. $ \mathbf{X} = [X_1 \; X_2 \;  ... \; X_n]^T $, $ \mathbf{1} = [1 \; 1 \; ... \; 1]^T $ 일 때 이를 선형대수적으로 다시 쓰면 다음과 같이 표현할 수 있다.


$$
\mathbf{X} \sim N(\mu \mathbf{1}, \sigma^2 \mathbf{I}_n) \tag{1}
$$

정규분포의 중요한 성질로, 정규분포를 선형변환해도 정규분포를 따른다는 것이 있다. 어떤 선형변환 $ \mathbf{P} $ 를 정규분포를 따르는 벡터 $ \mathbf{X} $ 에 적용하면

$$
\mathbf{PX} \sim N(\mu \mathbf{P1}, \sigma^2 \mathbf{PP}^T) \tag{2}
$$

<br>
<br>

### 2. 중심화 행렬 $ \mathbf{P} $

자, 이제 $ \mathbf{X} $ 를 이용해 $ S^2 $ 의 식을 만들 것이다. 결론부터 말하자면 $ (n-1) S^2 $ 은 $ \mathbf{X} $ 를 중심화한(성분의 평균이 0이 되도록한) 벡터의 노름의 제곱과 같다.

중심화 행렬 $ P $ 는 다음과 같이 쓸 수 있다. $ \mathbf{PX} $ 를 생각해보면 각 성분이 $ X_i - \bar{X} $ 가 됨을 쉽게 확인할 수 있다.

$$
\mathbf{P} = \mathbf{I}_n - \frac{1}{n}\mathbf{\mathbf{1}\mathbf{1}}^T \tag{3}
$$

이는 $ \frac{1}{\sqrt{n}} \mathbf{1} $ 을 법선벡터로 가지는 $ n-1 $ 차원 초평면으로의 정사영행렬이기도 하다. 정사영행렬은 idempotent 이고 대칭이기 때문에 좋은 성질들을 가지고 이는 이후 증명에 유용하게 쓰인다. (대칭행렬은 직교대각화가능하고, idempotent 이므로 eigenvalue로 0 또는 1만을 가진다)

<br>
<br>

### 3. 스펙트럼 정리를 이용한 $ \mathbf{P} $의 분석

$ \mathbf{P} $는 대칭행렬이므로, 직교대각화 가능하다. (스펙트럼 정리)

$$
\mathbf{P} = \mathbf{QDQ}^T, \quad

\mathbf{D} = 
\begin{bmatrix}
\mathbf{I}_{n-1} & 0 \\
0 & 0
\end{bmatrix} \tag{4}
$$

$ \mathbf{Q} $는 직교행렬, $ \mathbf{D} $는 대각행렬이다. $ rank(\mathbf{P}) = n-1 $이므로 $ \mathbf{D} $의 eigenvalue는 1이 중복도 n - 1, 0이 중복도 1을 가진다. 따라서 일반성을 잃지 않고 $ \mathbf{D} $의 대각선의 마지막 대각성분만 0이고 나머지는 1이라고 하자. 또한 1과 0에 대응되는 eigenspace를 생각해보면 $ E_0 = span\{{\mathbf{1}}\}, E_1 = {E_0}^{\perp} $ 임을 알 수 있다. 

<br>
<br>

### 4. $ \lVert \mathbf{PX} \rVert ^2 $ 의 계산

$$
\mathbf{(claim)} \quad (n-1) S^2 = \sigma^2 \sum_{i=1}^{n-1} Z_i^2
$$

위에서 $ (n-1) S^2 = \lVert \mathbf{PX} \rVert ^2 $ 임을 보였다.

$$
\lVert \mathbf{PX} \rVert ^2 = \lVert \mathbf{QDQ}^T\mathbf{X} \rVert ^2 = \lVert \mathbf{DQ}^T\mathbf{X} \rVert ^2 = \lVert \mathbf{DY} \rVert ^2
$$

은근슬쩍 $ \mathbf{Y} = \mathbf{Q}^T\mathbf{X} $ 치환을 했다. 이제 해야할 것이 좀 더 명확해졌다. $ \mathbf{DY} $ 에서는 $ \mathbf{Y} $ 의 앞의 n - 1개 성분만 살아남는다. 이들이 각각 정규분포를 따름을 보이면 된다!

수식 $ (2) $ 를 이용하면

$$
\mathbf{Q}^T\mathbf{X} \sim N(\mu \mathbf{Q}^T\mathbf{1}, \sigma^2 \mathbf{I}_n) \tag{5}
$$

$ \mathbf{Q} $는 직교행렬이므로, $ \mathbf{Q}^T\mathbf{Q} = \mathbf{I}_n$ 임은 바로 이용했다. $ (4) $에 의해, $ \mathbf{Q} = [\mathbf{q}_1 \; \mathbf{q}_2 \; ... \; \mathbf{q}_n]^T $ 일 때

$$
\mathbf{q}_1, \mathbf{q}_2, \cdots, \mathbf{q}_{n-1} \in span\{\mathbf{1}\}^{\perp}
$$


따라서 $ \mathbf{q}_1, \mathbf{q}_2, \cdots, \mathbf{q}_{n-1} $ 와 $ \mathbf{1} $ 의 내적은 0이다. 이를 이용해 평균 항을 계산하면


$$
\mathbf{Q}^T\mathbf{1} = 

\begin{bmatrix}
{\mathbf{q}_1}^T \\
{\mathbf{q}_2}^T \\
\vdots \\
{\mathbf{q}_{n-1}}^T \\
{\mathbf{q}_n}^T
\end{bmatrix}
\cdot \mathbf{1} = 

\begin{bmatrix}
{\mathbf{q}_1}^T\mathbf{1} \\
{\mathbf{q}_2}^T\mathbf{1} \\
\vdots \\
{\mathbf{q}_{n-1}}^T\mathbf{1} \\
{\mathbf{q}_n}^T\mathbf{1}
\end{bmatrix} = 

\begin{bmatrix}
0 \\
0 \\
\vdots \\
0 \\
{\mathbf{q}_n}^T\mathbf{1}
\end{bmatrix}
$$

$ \mathbf{DY} \sim N(\mathbf{0}, \sigma^2 \mathbf{I}_n) $ 을 증명했다! $ \lVert \mathbf{DY} \rVert ^2 = \sigma^2 \sum_{i=1}^{n-1} Z_i^2 $ 이 된다.

<br>
<br>

### 5. 조금 다른 방법

조금 관점이 다른 방법을 소개한다. 위에서는 그저 정직하게 계산을 했는데, 다음 정리를 이용하면 더 쉽게 생각할 수 있다.

***$ \mathbf{Z} \sim N(\mathbf{0}, \mathbf{I}_n) $ 일 때, $ \mathbf{Z} $ 를 $r$차원 부분공간에 사영한 벡터의 노름 제곱은 $ \chi^2(r) $을 따른다.***

위에서 했던 과정과 비슷하게 생각하면 쉽게 증명할 수 있다.

이제 $ \mathbf{X} \sim N(\mu \mathbf{1}, \sigma^2 \mathbf{I}_n) $ 일 때 $ \mathbf{X} $ 의 정규화를 생각하자. 

$$
\frac{\mathbf{X} - \mu \textbf{1}}{\sigma} = \mathbf{Z}
$$

$ \mathbf{X} $를 $ \mathbf{Z} $에 대해 정리해 $ \lVert \mathbf{PX} \rVert ^2 $ 에 대입 후 계산하면

$$
\lVert \mathbf{PX} \rVert ^2 = \lVert \mathbf{P(\mu \mathbf{1}+\sigma \mathbf{Z})} \rVert ^2 = \lVert \mathbf{\sigma PZ} \rVert ^2 = \sigma^2 \lVert \mathbf{PZ} \rVert ^2
$$

위의 명제에 의해 $ \lVert \mathbf{PZ} \rVert ^2 \sim \chi^2(n-1)$


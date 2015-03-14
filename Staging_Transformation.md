The first change of variables discussed is known as a _staging transformation_,

<a href='http://www.codecogs.com/eqnedit.php?latex=x_{1}=u_{1}'><img src='http://latex.codecogs.com/gif.latex?x_{1}=u_{1}%.png' title='x_{1}=u_{1}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=x_{i}=u_{i}@plus;\frac{i-1}{i}x_{i@plus;1}@plus;\frac{1}{i}u_{1}'><img src='http://latex.codecogs.com/gif.latex?x_{i}=u_{i}+\frac{i-1}{i}x_{i+1}+\frac{1}{i}u_{1}%.png' title='x_{i}=u_{i}+\frac{i-1}{i}x_{i+1}+\frac{1}{i}u_{1}' /></a> , where <a href='http://www.codecogs.com/eqnedit.php?latex=i=2,\ldots, P'><img src='http://latex.codecogs.com/gif.latex?i=2,\ldots, P%.png' title='i=2,\ldots, P' /></a>.

With this transformation the masses are no longer uniform. Instead they become,

<a href='http://www.codecogs.com/eqnedit.php?latex=m_{1}=0'><img src='http://latex.codecogs.com/png.latex?m_{1}=0%.png' title='m_{1}=0' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=m_{k}=\frac{k}{k-1}m,'><img src='http://latex.codecogs.com/png.latex?m_{k}=\frac{k}{k-1}m,%.png' title='m_{k}=\frac{k}{k-1}m,' /></a>   where    <a href='http://www.codecogs.com/eqnedit.php?latex=k=2,\ldots,P'><img src='http://latex.codecogs.com/png.latex?k=2,\ldots,P.%.png' title='k=2,\ldots,P' /></a>

so that the new masses become <a href='http://www.codecogs.com/eqnedit.php?latex=\tilde{m}_{1}=m'><img src='http://latex.codecogs.com/png.latex?\tilde{m}_{1}=m%.png' title='\tilde{m}_{1}=m' /></a> and <a href='http://www.codecogs.com/eqnedit.php?latex=\tilde{m}_{k}=\frac{k}{k-1}m'><img src='http://latex.codecogs.com/png.latex?\tilde{m}_{k}=\frac{k}{k-1}m%.png' title='\tilde{m}_{k}=\frac{k}{k-1}m' /></a> .

Furthermore, this transformation leads to the following equations of motion

<a href='http://www.codecogs.com/eqnedit.php?latex=\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}'><img src='http://latex.codecogs.com/gif.latex?\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}%.png' title='\dot{u}_{i}=\frac{p_{i}}{\tilde{m}_{i}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\dot{p}_{i}=-\frac{\partial U_{eff}^{stage}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}'><img src='http://latex.codecogs.com/gif.latex?\dot{p}_{i}=-\frac{\partial U_{eff}^{stage}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}%.png' title='\dot{p}_{i}=-\frac{\partial U_{eff}^{stage}}{\partial u_{i}}=-m_{i}\omega_{P}^{2}u_{i}-\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}' /></a>

Finally using the chain rule we arrive upon the following forces

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{i=1}^{P}\frac{\partial \phi}{\partial x_{i}}'><img src='http://latex.codecogs.com/gif.latex?\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{i=1}^{P}\frac{\partial \phi}{\partial x_{i}}%.png' title='\frac{1}{P}\frac{\partial \phi}{\partial u_{1}}=\frac{1}{P}\sum_{i=1}^{P}\frac{\partial \phi}{\partial x_{i}}' /></a>

<a href='http://www.codecogs.com/eqnedit.php?latex=\frac{\frac{1}{P}\partial \phi}{\partial u_{i}}=\frac{1}{P}[\frac{\partial \phi}{\partial x_{i}}@plus;\frac{i-2}{i-1}\frac{\partial \phi}{\partial x_{i-1}}]'><img src='http://latex.codecogs.com/gif.latex?\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}=\frac{1}{P}[\frac{\partial \phi}{\partial x_{i}}+\frac{i-2}{i-1}\frac{\partial \phi}{\partial x_{i-1}}]%.png' title='\frac{1}{P}\frac{\partial \phi}{\partial u_{i}}=\frac{1}{P}[\frac{\partial \phi}{\partial x_{i}}+\frac{i-2}{i-1}\frac{\partial \phi}{\partial x_{i-1}}]' /></a>.


---

**NEXT SECTION**: [Normal\_Mode\_Transformation](Normal_Mode_Transformation.md)

---

**SKIP TO:** [Introduction](Introduction.md), [Variable\_Transformations](Variable_Transformations.md), [Thermostats](Thermostats.md), [Results](Results.md), [Conclusions](Conclusions.md)
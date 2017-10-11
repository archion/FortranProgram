# lib 

##Dependence
- lapack95

##usage

### M_matrix
  - det_ratio_row(c,ci,A,iA,pb)
  - det_ratio_col(c,ci,A,iA,pb)
  - det_ratio_rowcol(c,ci,A,iA,iY,pb)
  - inv_update_row(c,ci,pb,A,iA)
  - inv_update_col(c,ci,pb,A,iA)
  - inv_update_rowcol(c,ci,iY,A,iA)

### M_lattice
```
latt:
输入
	a1: 基矢
	a2: 基矢
	T1: 边界矢量
	T2: 边界矢量
	bdc: 边界条件
	layer: 层数
	sub: 子格的空间坐标
输出
	Ns: 总格点数
	i2r: 格点指标对应的坐标
	neb: 第几个格点的近邻
		nb: 第几近邻
			neb: 近邻格点的格点指标
			bond: 近邻键的键指标
			bdc: 跳跃附带的边界相位
			dr: 近邻键的方向向量
	bond: 第几近邻键
		bd: 键指标
			i: 键的两个格点指标
			r: 键中心坐标
			dr: 键的方向向量
			bdc: 键的边界相位
方法
	gen_latt(): 生成 latt%i2r
	gen_neb(l): 生成 latt%neb 到 l 近邻
	gen_bond(l): 生成 latt%bond  到 l 近邻 
	gen_brizon(brizon): 生成 brizon
	gen_origin_brizon(a1,a2,o_brizon): 以当前brizon生成正格矢为a1，a2的布里渊区o_brizon

brizon:
输入
	n1: k点数目
	n2: k点数目
输出
	b1: 倒格矢
	b2: 倒格矢
	k: k点的坐标
	T: 布里渊区的端点坐标

check_lattice(ut): 输出晶格到文件ut
```

### M_lattice_final
```
latt:
输入
	is_all: 是否生成整个晶格（否则只生成super cell，默认.false.）
	a1: 基矢
	a2: 基矢
	c1: super cell基矢
	c2: super cell基矢
	T1: 边界矢量
	T2: 边界矢量
	bdc: 边界条件
	layer: 层数（默认1，多层还没实现）
	rsb: 子格的空间坐标
	i2isb(i,sb): 第i个格点之前（包括当前格点）有多少个sb子格 
	isb2i(i,sb): 第sb子格的第i个格点对应的格点指标
输出
	sb: 子格数
	Ns: 总格点数
	Nc: super cell个数
	Ni: super cell内的格点数
	nb: 近邻指标
		bd: 键指标
			sb: 两个格点的子格指标
			sbc: 两个格点的super cell指标 
			i: 两个格点的格点指标（若不生成整个晶格，则格点指标和super cell指标相同）
			r: 键的中心坐标
			dr: 键的方向向量（由i(1)指向i(2)）
			bdc: 键上带的边界条件 
		st: 格点指标
			j: 近邻的格点指标
			bd: 近邻的键指标
			dir: 近邻的方向是否与键方向一致
方法
	gen_latt(l): 生成包含l近邻的latt
	gen_brizon(brizon): 生成brizon

brizon:
输出
	nk: 独立k点数
	nq: 非独立的k点对应的平移动量数
	a1: 格子的倒格矢
	a2: 格子的倒格矢
	Ta: 格子对应的第一布里渊区端点坐标
	c1: super cell的倒格矢
	c2: super cell的倒格矢
	Tc: super cell对应的第一布里渊区端点坐标
	k0: 由边界条件导致的k点的移动
	k: 全部的k点坐标
	q: 非独立的k点对应的平移动量 

check_lattice(ut): 输出晶格到文件ut
```

### M_hamilton
```
输入
	spin: 自旋数，默认为2
	is_sc: 是否有超导，默认为.true.
输出
	var(:): 不同类型的变分参数，其中var(:0)为不变分，var(1:)做变分
		val(:)：变分参数的值
		v2i(:)：对应变分参数的格点或键指标
			i(:)：格点或键指标
		i2v(:)：对应格点或键所对应的变分参数指标
		bd_sg(:)：对应格点或键指标所带的因子
		n：该类型变分参数的个数
		V：该类型变分参数所带的相互作用
		Vn：相互作用是第几近邻
		nb：该类型变分参数是第几近邻
		sg：变分参数的类型
	方法
		Hamilton：生成平均场Hamiltonian
		dHamilton：生成平均场Hamiltonian对变分参数的导数
		get：按val向量设置变分参数的值
		put：生成保护变分参数值的向量val

gen_var(sg,nb,V,val,Vn):
var_shrink()
Green(gm,k,omg)
LDOS(ut,gm,i,omg,m)
EDC(ut,gm,k,omg,m)
DOS(ut,gm,omg,m,peak)
fermis(ut,gm,k,omg)
energy(ut,k)
band_e(ut,gm,ki,kf,n,omg,m)
band(ut,ki,kf,n)
```

### M_hamilton_final
```
输入
	spin: 自旋数，默认为2
	is_nambu: 是否写到Nambu表象，默认为.true.
	nf: 粒子数
	Tk: 温度
输出
	Uqi: 从supercell的格点指标到q的Fourier变换（Uqi*c(k,i)=c(k+q)），顺序是：子格，q，自旋
	Uik: 从独立的k点到实空间的Fourier变换（Uik*c(k,i)=c(i)）
	var(:): 不同类型的变分参数，其中var(:0)为不变分，var(1:)做变分
		val(:)：变分参数的值
		v2i(:)：对应变分参数的格点或键指标
			i(:)：格点或键指标
		i2v(:)：对应格点或键所对应的变分参数指标
		bd_sg(:)：对应格点或键指标所带的因子
		n：该类型变分参数的个数
		V：该类型变分参数所带的相互作用
		Vn：相互作用是第几近邻
		nb：该类型变分参数是第几近邻
		tp：变分参数的类型
	方法
		Hamilton：生成平均场Hamiltonian
		dHamilton：生成平均场Hamiltonian对变分参数的导数
		get：按val向量设置变分参数的值
		put：生成保护变分参数值的向量val

gen_var(tp,nb,V,val,Vn):
var_init()
Green(gm,k,omg)
LDOS(ut,gm,i,omg,m)
EDC(ut,gm,k,omg,m)
DOS(ut,gm,omg,m,peak)
fermis(ut,gm,k,omg)
energy(ut,k)
band_e(ut,gm,ki,kf,n,omg,m)
band(ut,ki,kf,n)
```

### M_hamilton_final_M
```
输入
	nf: 粒子数
	Tk: 温度
输出
	t_ham:
		fer_bos: 默认为费米子-1，玻色子为1
		is_nambu: 是否为Nambu表象
		mask(:,sb): 该态排序指标在子格上是否有态
		i2tp(:): 态排序的指标到态指标的映射（厄米共轭对应于不同的态指标，用±号表示）
		tp2i(:): 态的指标到态排序指标的映射
		ptp(:): 态排序指标之前（包括当前指标）有多少个supercell格点指标（乘以supercell个数可以得到格点指标）
		i2H(:,tp): supercell格点指标和态指标到平均场Hamilton矩阵指标的映射
		s2H(:,tp): 格点指标和态指标到平均场Hamilton矩阵指标的映射
		H2i(:,2): 平均场Hamilton矩阵指标到supercell格点指标和态指标的映射
		H2s(:,2): 平均场Hamilton矩阵指标到格点指标和态指标的映射
		Hi: 在supercell表象下平均场Hamilton矩阵的大小
		Hs: 在格点表象下平均场Hamilton矩阵的大小
		Uqi(q,i): 从supercell的表象的Hamiltonian指标到q的Fourier变换（Uqi*c(k,i)=c(k+q)），顺序是：子格，q，自旋
		Uik(i,k): 从supercell的格点表象的Hamiltonian指标到格点表象的Hamiltonian指标的Fourier变换（Uik*c(k,i)=c(i)）
		rg(2): var指标的上下限
		var(:): 不同类型的变分参数，其中var(:0)为不变分，var(1:)做变分
			nb: 变分参数定义在第几近邻
			c(:,:): 该变分参数耦合的产生湮灭算符，第二指标表示不同项 
			sg: 不同相之间的相对系数
			cg: 不同相是否需要对bd取复共轭
			V: 相互作用大小，不设置则表示不做变分
			val(:): 变分参数的值
			bd(:)：对应格点或键指标所带的因子
			bd2v(:)：对应格点或键所对应的变分参数指标
			n：该类型变分参数的个数
			label: 该类型变分参数的字符标志
		方法
			add: 添加Hamilton项
			init: 初始化Hamiltion量
			Hamilton：生成平均场Hamiltonian
			dHamilton：生成平均场Hamiltonian对变分参数的导数
			get：按val向量设置变分参数的值
			put：生成保护变分参数值的向量val
			Green(gm,k,omg)
			LDOS(ut,gm,i,omg,m)
			EDC(ut,gm,k,omg,m)
			DOS(ut,gm,omg,m,peak)
			fermis(ut,gm,k,omg)
			energy(ut,k)
			band_e(ut,gm,ki,kf,n,omg,m)
			band(ut,ki,kf,n)

c:
	l: "i" 
	sb: 子格指标
	tp: 态指标
```


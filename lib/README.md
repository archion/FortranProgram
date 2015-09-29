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

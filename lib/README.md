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
	gen_brizon(): 生成 brizon

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


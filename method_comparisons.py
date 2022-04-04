import numpy as np
import ssDMEWMAC.functions as spc

# a sample data frame from Hawkins & Maboudou-Tchao (2007)
X = list((np.array([-1.347, 25.600, 12.400, 0.256, 58.600]),
          np.array([-1.386, 24.800, 13.300, 0.285, 58.100]),
          np.array([-0.478, 24.900, 12.700, 0.248, 57.800]),
          np.array([-1.772, 25.300, 12.200, 0.265, 58.200]),
          np.array([-1.309, 25.100, 12.300, 0.236, 58.800]),
          np.array([-0.580, 24.500, 12.200, 0.252, 56.700]),
          np.array([-0.892, 26.000, 12.200, 0.270, 58.600]),
          np.array([-1.470, 24.400, 12.700, 0.248, 57.900]),
          np.array([-0.673, 25.600, 12.100, 0.218, 58.100]),
          np.array([-0.315, 25.400, 12.400, 0.212, 57.900])
          ))
factor_matrix = np.zeros([6,6], dtype=np.float)
n = 1
m_prior = np.zeros(5)
d_prior = np.zeros(5)

for x in X:
    factor_matrix, r, u = spc.u_transform(factor_matrix=factor_matrix, obs_vector=x, n=n, p=5)

    if n == 7:
        m = spc.m_n(smoothing=0.1,u_n=u,m_prior=m_prior)
        t2, mn2 = spc.t_n(smoothing=0.1,m=m,n=n,p=5)
        d = spc.d_n(smoothing=0.1, m=m, d_prior=d_prior)
        dt2, dn2 = spc.dt_n(smoothing=0.1, d=d, n=n, p=5)
    elif n > 7:
        m_prior = m
        m =  spc.m_n(smoothing=0.1,u_n=u,m_prior=m_prior)
        t2, mn2 = spc.t_n(smoothing=0.1, m=m, n=n, p=5)
        d_prior = d
        d = spc.d_n(smoothing=0.1, m=m, d_prior=d_prior)
        dt2, dn2 = spc.dt_n(smoothing=0.1, d=d, n=n, p=5)
    else:
        m = np.repeat(np.nan, 5, axis=0)
        t2 = np.nan
        mn2 = np.nan
        d = np.repeat(np.nan, 5, axis=0)
        dt2 = np.nan
        dn2 = np.nan
    n += 1
    #print(factor_matrix)
    print(r)
    print(u)
    print(m)
    #print(t2)
    print(mn2)
    #print(d)
    #print(dt2)
    #print(dn2)


# now let us make sure we duplicated the DMEWMA code properly
X = list((
    np.array([-1.19, 0.59]),
    np.array([0.12, 0.90]),
    np.array([-1.69, 0.40]),
    np.array([0.30, 0.46]),
    np.array([0.89, -0.75]),
    np.array([0.82, 0.98]),
    np.array([-0.30, 2.28]),
    np.array([0.63, 1.75]),
    np.array([1.56, 1.58]),
    np.array([1.46, 3.05])
))

n = 1
y_prior = np.zeros(2)
z_prior = np.zeros(2)

for x in X:
    if n == 1:
        y = spc.y_vec(smoothing=0.1,x_n=x,y_prior=y_prior)
        z = spc.z_vec(smoothing=0.1,y_n=y,z_prior=z_prior)
        tm2n = spc.mewma(smoothing=0.1,n=n,yvec=y,sigma_0=np.array([[1,0.5],[0.5,1]]))
        td2n = spc.dmewma(smoothing=0.1,n=n,zvec=z,sigma_0=np.array([[1,0.5],[0.5,1]]))
    if n > 1:
        y_prior = y
        z_prior = z
        y = spc.y_vec(smoothing=0.1, x_n=x, y_prior=y_prior)
        z = spc.z_vec(smoothing=0.1, y_n=y, z_prior=z_prior)
        tm2n = spc.mewma(smoothing=0.1, n=n, yvec=y, sigma_0=np.array([[1, 0.5], [0.5, 1]]))
        td2n = spc.dmewma(smoothing=0.1, n=n, zvec=z, sigma_0=np.array([[1, 0.5], [0.5, 1]]))
    n += 1
    #print(z)
    #print(td2n)
    print(y)
    print(tm2n)


#arl = simulate('MEWMA',2,0.1,8.66,100000)
#print(arl.mean())
# 196.0808

#arl2 = simulate('ssMEWMA',2,0.1,8.786,100000)
#print(arl2.mean())
#204.6277



h_mewma_2_1 = spc.find_h('MEWMA', 2, 0.1, 500, 10.827, 0.001, 1000)
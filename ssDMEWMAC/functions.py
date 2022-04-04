"""
This is part of the python program that will be used to run the derivation and
the performance of the ssDMEWMAC control chart. It is a conversion of the work
of the ssMEWMAC Fortran 95 code into something more user friendly to run.

The key principles of the routine are to perform the Cholesky decomposition of
the covariance matrix, pass it to the converted ssMEWMA routine provided in
Hawkins and Maboudou-Tchao (2007). Within this routine, the factored matrix of
the covariance matrix is used to calculate recursive residuals of which are
then transformed into multivariate standard normal vectors. The standardized
vectors of recursive residuals are the primary values passed to the ssMEWMA &
ssDMEWMAC test statistics.
"""

import numpy as np
from scipy.stats import t
from scipy.stats import norm
from math import copysign
from sklearn.datasets import make_sparse_spd_matrix
from google.cloud import storage


# the critical u transformation function
def u_transform(factor_matrix, obs_vector, n, p, epsilon=1e-20):
    """
    :param factor_matrix: cholesky inverse root matrix of the covariance matrix
    :param obs_vector: p-dimensional observation vector
    :param n: the sample n
    :param p: the dimensionality of the process
    :param epsilon: the tolerance value of the function
    :return: updated factor_matrix, recursive residual vector, and u-vector
    """

    # TODO: Checks to make sure factor_matrix and obs_vector match p dimensions

    # if n == 1, initialize the factor matrix to 0 matrix
    if n == 1:
        factor_matrix = np.zeros([p + 1, p + 1], dtype=np.float)
    else:
        factor_matrix = factor_matrix

    # initialize a zero vector for R and a missing values vector for U
    r = np.zeros(p, dtype=np.float)
    u = np.repeat(np.nan, p)

    # initialize a copy of the diagonal elements of the factor_matrix
    old = np.repeat(np.nan, p + 1)

    # Now to translate the existing Fortran code to Python...
    # step one, take the observation vector and prepend a 1 to it
    work = np.insert(obs_vector, 0, 1.0, axis=0)

    # Now we are going to create the factor_matrix numerically
    for i in range(p + 1):
        d1 = factor_matrix[i, i]
        old[i] = d1
        d2 = work[i]
        if i > 0:
            r[i - 1] = d2
        d0 = max(np.abs(d1), np.abs(d2))
        if d0 != 0:
            d1 = d1 / d0
            d2 = d2 / d0
            if np.abs(d2) >= epsilon:
                d = np.sqrt((d1 ** 2 + d2 ** 2))
                factor_matrix[i, i] = d * d0
                c0 = d1 / d
                s0 = d2 / d
                for j in range((i + 1), p + 1):
                    d1 = factor_matrix[i, j]
                    d2 = work[j]
                    factor_matrix[i, j] = d1 * c0 + d2 * s0
                    work[j] = -d1 * s0 + d2 * c0
    for i in range(1, p + 1):
        ndf = n - i - 1
        df = ndf
        if ndf > 0 and old[i] > epsilon:
            tee = r[i - 1] * np.sqrt(df) / old[i]
            if ndf < 13:
                area = t.cdf(tee, ndf)
                # this gets the inverse normal value of the area calculated
                u[i - 1] = norm.ppf(area)
            else:
                sign = copysign(1, tee)
                u[i - 1] = sign * (((1.0 - 2.0 / (8.0 * df + 3.0)) * \
                                    np.sqrt(df * np.log(1 + tee ** 2 / df))))

    return factor_matrix, r, u


# function for calculating the M vector
def m_n(smoothing, u_n, m_prior):
    """
    :param smoothing: value of the smoothing constant lambda
    :param u_n: the u_n vector from the u_transform function
    :param m_prior: the previous output of this function or 0 vector size p
    :return: m_n vector of size p (will occur at p + 2 sample)
    """

    m = smoothing * u_n + (1 - smoothing) * m_prior

    return m


# function for calculating the plotting statistic
def t_n(smoothing, m, n, p):
    """
    :param smoothing: value of the smoothing constant lambda
    :param m: the m vector at sample n
    :param n: the sample n
    :param p: the dimensionality of the process
    :return: the t^2 statistic for plotting
    """

    # need to calculate Sigma_m
    Sigma_m = (smoothing / (2 - smoothing)) * \
              (1 - (1 - smoothing) ** (2 * (n - p - 1))) * \
              np.identity(p)

    # get the inverse matrix
    Sigma_m_inv = np.linalg.inv(Sigma_m)

    # now calculate the t^2 statistic
    t = m.dot(Sigma_m_inv).dot(m)

    # and mn2 for plotting convenience
    mn2 = m.dot(m)

    return t, mn2


# function for plotting the ssMEWMA bound at sample n
def ssmewma_bound(smoothing, n, p, h):
    """
    :param smoothing: value of the smoothing constant lambda
    :param n: the sample n
    :param p: the dimensionality of the process
    :param h: the value h to control the ARL0
    :return: the value of the confidence bound at sample n
    """
    exponent_term = (1 - (1 - smoothing) ** (2 * (n - p - 1)))

    c = ((smoothing * exponent_term) / (2 - smoothing)) * h

    return c


# function for calculating the ssMEWMC control statistic
def s_n(smoothing, u_n, s_prior):
    """
    :param smoothing: value of the smoothing constant lambda
    :param u_n: the u-vector at time n
    :param s_prior: the prior value from the s_n function
    :return: the current s matrix at time n
    """
    s = (1 - smoothing) * s_prior + smoothing * (np.outer(u_n, u_n.T))

    return s


# function for the plotting statistic of the ssMEWMC
def c_n(s, p):
    """
    :param s: the s_n out put at time n
    :param p: the dimensionality of the process
    :return: the plotting statistic
    """
    # TODO: Need to fix this for when the determinate is 0
    c = s.trace() - np.log(np.linalg.det(s)) - p

    return c


# define the relevant fucntions for the ssDMEWMAC charts
# start with the ssDMEWMA set of functions
def d_n(smoothing, m, d_prior):
    """
    :param smoothing: value of the smoothing constant lambda
    :param m: value of ssMEWMA m vector at sample n
    :param d_prior: the prior value of this function
    :return: d vector at sample n
    """
    d = smoothing * m + (1 - smoothing) * d_prior

    return d


# function to get the ssDMEWMA plotting statistics
def dt_n(smoothing, d, n, p):
    """
    :param smoothing: value of the smoothing constant lambda
    :param d: the d vector at sample n
    :param n: the sample n
    :param p: the dimensionality of the process
    :return: the t^2 statistic for plotting
    """

    # need to calculate Sigma_m
    quad = (smoothing ** 4)
    sq = (1 - smoothing) ** 2
    exp1 = ((n + 1) ** 2) * ((1 - smoothing) ** (2 * n))
    exp2 = ((2 * n ** 2) + (2 * n) - 1) * ((1 - smoothing) ** (2 * n + 2))
    exp3 = ((n ** 2) * ((1 - smoothing) ** (2 * n + 4)))
    denom = (1 - (1 - smoothing) ** 2) ** 3

    Sigma_d = ((quad * (1 + sq - exp1 + exp2 - exp3)) / denom) * \
              np.identity(p)

    # get the inverse matrix
    Sigma_d_inv = np.linalg.inv(Sigma_d)

    # now calculate the t^2 statistic
    t = d.dot(Sigma_d_inv).dot(d)

    # and mn2 for plotting convenience
    dn2 = d.dot(d)

    return t, dn2


# function to calculate the ssDMEWMA bound at sample n
def dssmewma_bound(smoothing, n, h):
    """
    :param smoothing: value of the smoothing constant lambda
    :param n: the sample n
    :param h: the value h to control the ARL0
    :return: the value of the confidence bound at sample n
    """
    quad = (smoothing ** 4)
    sq = (1 - smoothing) ** 2
    exp1 = ((n + 1) ** 2) * ((1 - smoothing) ** (2 * n))
    exp2 = ((2 * n ** 2) + (2 * n) - 1) * ((1 - smoothing) ** (2 * n + 2))
    exp3 = ((n ** 2) * ((1 - smoothing) ** (2 * n + 4)))
    denom = (1 - (1 - smoothing) ** 2) ** 3

    c = ((quad * (1 + sq - exp1 + exp2 - exp3)) / denom) * h

    return c


# function to calculate the c_n metric
def v_n(smoothing, s, v_prior):
    """
    :param smoothing: value of the smoothing constant lambda
    :param s: the s matrix at sample n
    :param v_prior: the prior value from this function
    :return: the v matrix at sample n
    """
    v = smoothing * s + (1 - smoothing) * v_prior

    return v


# the control statistic for the ssDMEWMC
def cv_n(v, p):
    """
    :param v: the v matrix at sample n
    :param p: the dimensionality of the process
    :return: the value of the control statistic at sample n
    """
    # TODO: Need to fix this for when the determinate is 0
    c = v.trace() - np.log(np.linalg.det(v)) - p

    return c


# need a set of functions to replicate the original (D)MEWMA
def y_vec(smoothing, x_n, y_prior):
    """
    :param smoothing: value of the smoothing constant lambda
    :param x_n: the value of the x vector at sample n
    :param y_prior: the prior value from this function
    :return: the value of the y vector at sample n 
    """
    return smoothing * x_n + (1 - smoothing) * y_prior


def z_vec(smoothing, y_n, z_prior):
    """
    :param smoothing: value of the smoothing constant lambda
    :param y_n: the value of the y vector at sample n
    :param z_prior: the prior value from this function
    :return: the value of the z vector at sample n 
    """
    return smoothing * y_n + (1 - smoothing) * z_prior


def mewma(smoothing, n, yvec, sigma_0):
    """
    :param smoothing: value of the smoothing constant lambda
    :param n: the sample n
    :param y_vec: the y_vector at sample n
    :param sigma_0: the value of the initial covariance matrix
    :return: the value of the MEWMA algorithm at sample n
    """

    sigma_y = ((smoothing * (1 - (1 - smoothing) ** (2 * n))) / (2 - smoothing)) * sigma_0

    t2_n = yvec.dot(np.linalg.inv(sigma_y)).dot(yvec)

    return t2_n


def dmewma(smoothing, n, zvec, sigma_0):
    """
    :param smoothing: value of the smoothing constant lambda
    :param n: the sample n
    :param z_vec: the value of the z vector at sample n
    :param sigma_0: the value of the initial covariance matrix
    :return: the value of the DMEWMA algorithm at sample n
    """

    quad = (smoothing ** 4)
    sq = (1 - smoothing) ** 2
    exp1 = ((n + 1) ** 2) * ((1 - smoothing) ** (2 * n))
    exp2 = ((2 * n ** 2) + (2 * n) - 1) * ((1 - smoothing) ** (2 * n + 2))
    exp3 = ((n ** 2) * ((1 - smoothing) ** (2 * n + 4)))
    denom = (1 - (1 - smoothing) ** 2) ** 3

    sigma_z = ((quad * (1 + sq - exp1 + exp2 - exp3)) / denom) * sigma_0

    tz2_n = zvec.dot(np.linalg.inv(sigma_z)).dot(zvec)

    return tz2_n


# function for creating random mean and covariance matrices size p
def sampling_distribution(p, seed, alpha=0.5, norm_diag=True,
                          smallest_coef=0.1, largest_coef=2):
    cov_mat = make_sparse_spd_matrix(dim=p,
                                     alpha=alpha,
                                     norm_diag=norm_diag,
                                     smallest_coef=smallest_coef,
                                     largest_coef=largest_coef,
                                     random_state=seed)

    mean_vec = np.zeros(p, dtype=float)

    """
    for i in range(p):
        subseed = seed*i
        random_percent = np.random.uniform(0.01,0.20,1)
        variance = cov_mat[i,i]
        sd = np.sqrt(variance)
        mean_vec[i] = sd / random_percent
    """

    return mean_vec, cov_mat


# function for running the simulation
def simulate(method, dimension, smoothing, h, iterations, out_of_control=False, delta=0.5, change_after=1):
    """
    :param mean_vector: the multivariate mean vector to sample from
    :param covariance_matrix: the covariance matrix to sample from
    :param method: either "ssMEWMAC" or "ssDMEWMAC" method
    :param h: the h-value to test for o.o.c signaling
    :param iterations: the number of iterations to perform
    :param out_of_control: Boolean indicator if an out of control simulation is taking place; default false
    :param delta: The magnitude of the out of control signal to simulate
    :param change_after: The number of samples to take before inducing the out of control signal
    :return: the in control average run length for this arrangement
    """

    # deduce p
    p = dimension

    # hold the arl's
    arl = np.repeat(np.nan, iterations, axis=0)

    for i in range(iterations):
        # print(i)
        j = int((i + 1))

        # initialize relevant variables
        factor_matrix = np.zeros([p + 1, p + 1], dtype=np.float)
        n = 1
        m_prior = np.zeros(p)
        d_prior = np.zeros(p)
        s_prior = np.identity(p)
        v_prior = np.identity(p)
        y_prior = np.zeros(p)
        z_prior = np.zeros(p)

        # hold the vector arrays in a list for archiving
        X = list()
        R = list()
        U = list()
        M = list()
        D = list()
        mewmas = list()
        mns = list()
        dmewmas = list()
        dns = list()
        cns = list()
        cvns = list()

        mn2 = 0
        dn2 = 0
        cn = 0
        cvn = 0

        np.random.seed(j)

        if method == "ssMEWMC" or method == "ssDMEWMC":
            norm = False
        else:
            norm = True

        mean_vector, covariance_matrix = sampling_distribution(p=p, seed=j,
                                                               alpha=np.random.uniform(0.25, 0.75, 1),
                                                               norm_diag=norm, smallest_coef=0.1,
                                                               largest_coef=1)

        obs = np.random.multivariate_normal(mean_vector, covariance_matrix, 20000)

        if out_of_control != False:
            t = p + change_after

            if method != "ssMEWMC" and method != "ssDMEWMC":
                # need to solve delta = (mu' Sigma^-1 mu)^(1/2) for mu
                d = delta

                # for proof of concept, and without loss of generality, add d to every row after time t to simulate o.o.c
                mean_vector[0] = mean_vector[0] + d
                obs_ooc = np.random.multivariate_normal(mean_vector, covariance_matrix, 20000 - t + 1)

                obs = np.concatenate((obs[0:t - 1, :], obs_ooc))
            if method == "ssMEWMC" or method == "ssDMEWMC":
                covariance_matrix[0, 0] = delta * covariance_matrix[0, 0]

                obs_ooc = np.random.multivariate_normal(mean_vector, covariance_matrix, 20000 - t + 1)

                obs = np.concatenate((obs[0:t - 1, :], obs_ooc))

        if method == "ssMEWMA":
            while n <= 20000:
                # print(n)
                ob = obs[n - 1]

                # print(ob)
                X.append(ob)

                factor_matrix, r, u = u_transform(factor_matrix=factor_matrix, obs_vector=ob, n=n, p=p, epsilon=1e-20)
                R.append(r)
                U.append(u)

                if n == p + 2:
                    m = m_n(smoothing=smoothing, u_n=u, m_prior=m_prior)
                    t2, mn2 = t_n(smoothing=smoothing, m=m, n=n, p=p)
                    M.append(m)
                    mns.append(mn2)
                elif n > p + 2:
                    m_prior = m
                    m = m_n(smoothing=smoothing, u_n=u, m_prior=m_prior)
                    t2, mn2 = t_n(smoothing=smoothing, m=m, n=n, p=p)
                    M.append(m)
                    mns.append(mn2)
                else:
                    m = np.repeat(np.nan, p, axis=0)
                    t2 = 0
                    mn2 = 0
                    M.append(m)
                    mns.append(mn2)

                if n > p + 2:
                    break_val = ssmewma_bound(smoothing, n, p, h)

                if n > p + 2 and mn2 > break_val:
                    arl[i] = n
                    break

                # prevent nulls
                if n == 20000:
                    arl[i] = n
                    break

                # print(mn2)

                j += 1
                n += 1

        elif method == "ssMEWMC":
            while n <= 20000:
                # print(n)
                ob = obs[n - 1]

                # print(ob)
                X.append(ob)

                factor_matrix, r, u = u_transform(factor_matrix=factor_matrix, obs_vector=ob, n=n, p=p, epsilon=1e-20)
                R.append(r)
                U.append(u)

                if n == p + 2:
                    s = s_n(smoothing=smoothing, u_n=u, s_prior=s_prior)
                    cn = c_n(s, p)
                    cns.append(cn)
                elif n > p + 2:
                    s_prior = s
                    s = s_n(smoothing=smoothing, u_n=u, s_prior=s_prior)
                    cn = c_n(s, p)
                    cns.append(cn)
                else:
                    s = np.identity(p)
                    cn = 0
                    cns.append(cn)

                if n > p + 2:
                    break_val = h

                if n > p + 2 and cn > break_val:
                    arl[i] = n
                    break

                # prevent nulls
                if n == 20000:
                    arl[i] = n
                    break

                # print(cn)

                j += 1
                n += 1

        elif method == "ssDMEWMA":
            while n <= 20000:
                # print(n)
                ob = obs[n - 1]

                # print(ob)
                X.append(ob)

                factor_matrix, r, u = u_transform(factor_matrix=factor_matrix, obs_vector=ob, n=n, p=p, epsilon=1e-20)
                R.append(r)
                U.append(u)

                if n == p + 2:
                    m = m_n(smoothing=smoothing, u_n=u, m_prior=m_prior)
                    t2, mn2 = t_n(smoothing=smoothing, m=m, n=n, p=p)
                    M.append(m)
                    mns.append(mn2)
                    d = d_n(smoothing=smoothing, m=m, d_prior=d_prior)
                    D.append(d)
                    dt2, dn2 = dt_n(smoothing, d, n, p)
                    dns.append(dn2)
                elif n > p + 2:
                    m_prior = m
                    m = m_n(smoothing=smoothing, u_n=u, m_prior=m_prior)
                    t2, mn2 = t_n(smoothing=smoothing, m=m, n=n, p=p)
                    M.append(m)
                    mns.append(mn2)
                    d_prior = d
                    d = d_n(smoothing=smoothing, m=m, d_prior=d_prior)
                    D.append(d)
                    dt2, dn2 = dt_n(smoothing, d, n, p)
                    dns.append(dn2)
                else:
                    m = np.repeat(np.nan, p, axis=0)
                    t2 = 0
                    mn2 = 0
                    M.append(m)
                    mns.append(mn2)
                    d = np.repeat(np.nan, p, axis=0)
                    dt2 = 0
                    dn2 = 0
                    D.append(d)
                    dns.append(dn2)

                if n > p + 2:
                    break_val = dssmewma_bound(smoothing, n, h)

                if n > p + 2 and dn2 > break_val:
                    arl[i] = n
                    break

                # prevent nulls
                if n == 20000:
                    arl[i] = n
                    break

                # print(dn2)

                j += 1
                n += 1

        elif method == "ssDMEWMC":
            while n <= 20000:
                # print(n)
                ob = obs[n - 1]

                # print(ob)
                X.append(ob)

                factor_matrix, r, u = u_transform(factor_matrix=factor_matrix, obs_vector=ob, n=n, p=p, epsilon=1e-20)
                R.append(r)
                U.append(u)

                if n == p + 2:
                    s = s_n(smoothing=smoothing, u_n=u, s_prior=s_prior)
                    cn = c_n(s, p)
                    cns.append(cn)
                    v = v_n(smoothing=smoothing, s=s, v_prior=v_prior)
                    cvn = cv_n(v=v, p=p)
                    cvns.append(cvn)
                elif n > p + 2:
                    s_prior = s
                    s = s_n(smoothing=smoothing, u_n=u, s_prior=s_prior)
                    cn = c_n(s, p)
                    cns.append(cn)
                    v_prior = v
                    v = v_n(smoothing=smoothing, s=s, v_prior=v_prior)
                    cvn = cv_n(v=v, p=p)
                    cvns.append(cvn)
                else:
                    s = np.identity(p)
                    cn = 0
                    cns.append(cn)
                    v = np.identity(p)
                    cvn = 0
                    cvns.append(cvn)

                if n > p + 2:
                    break_val = h

                if n > p + 2 and cvn > break_val:
                    arl[i] = n
                    break

                # prevent nulls
                if n == 20000:
                    arl[i] = n
                    break

                # print(cn)

                j += 1
                n += 1

        elif method == 'DMEWMA':
            while n <= 20000:
                # print(n)
                ob = obs[n - 1]

                # print(ob)
                X.append(ob)

                if n == 1:
                    y = y_vec(smoothing=smoothing, x_n=ob, y_prior=y_prior)
                    z = z_vec(smoothing=smoothing, y_n=y, z_prior=z_prior)
                    td2n = dmewma(smoothing=smoothing, n=n, zvec=z, sigma_0=covariance_matrix)
                if n > 1:
                    y_prior = y
                    z_prior = z
                    y = y_vec(smoothing=smoothing, x_n=ob, y_prior=y_prior)
                    z = z_vec(smoothing=smoothing, y_n=y, z_prior=z_prior)
                    td2n = dmewma(smoothing=smoothing, n=n, zvec=z, sigma_0=covariance_matrix)

                dmewmas.append(td2n)

                if n > p + 2:
                    break_val = h

                if n > p + 2 and td2n > break_val:
                    arl[i] = n
                    break

                # prevent nulls
                if n == 20000:
                    arl[i] = n
                    break

                # print(td2n)

                j += 1
                n += 1

        elif method == 'MEWMA':
            while n <= 20000:
                # print(n)
                ob = obs[n - 1]

                # print(ob)
                X.append(ob)

                if n == 1:
                    y = y_vec(smoothing=smoothing, x_n=ob, y_prior=y_prior)
                    tm2n = mewma(smoothing=smoothing, n=n, yvec=y, sigma_0=covariance_matrix)
                if n > 1:
                    y_prior = y
                    y = y_vec(smoothing=smoothing, x_n=ob, y_prior=y_prior)
                    tm2n = mewma(smoothing=smoothing, n=n, yvec=y, sigma_0=covariance_matrix)

                mewmas.append(tm2n)

                if n > p + 2:
                    break_val = h

                if n > p + 2 and tm2n > break_val:
                    arl[i] = n
                    break

                # prevent nulls
                if n == 20000:
                    arl[i] = n
                    break
                # print(tm2n)

                j += 1
                n += 1

        else:
            print("Please specify either 'MEWMA', 'DMEWMA', 'ssMEWMA/C' or 'ssDMEWMA/C' in order to proceed.")

    return arl


# function to upload to google cloud bucket
def upload_blob(bucket_name, source_file_name, destination_blob_name):
    """Uploads a file to the bucket."""
    # bucket_name = "your-bucket-name"
    # source_file_name = "local/path/to/file"
    # destination_blob_name = "storage-object-name"

    storage_client = storage.Client(project="fresch-analytics-283821")
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(destination_blob_name)

    with open(source_file_name) as f:
        first_line = f.readline()
        blob.upload_from_string(first_line)

    print(
        "File {} uploaded to {}.".format(
            source_file_name, destination_blob_name
        )
    )


# need a function then, that increments h_init, up or down to get the h_final to achieve the desired ARL
def find_h(method, dimension, smoothing, target_arl, h_init, tolerance, iterations):
    h_list = list()

    c = 0
    delta = target_arl * tolerance

    diff = delta + 0.01

    print('Tolerance = %s' % delta)

    while diff > delta:
        print(h_init)

        arl = simulate(method=method, dimension=dimension, smoothing=smoothing, h=h_init, iterations=iterations)
        arl0 = np.nanmean(arl)
        print(arl0)

        if c == 0:
            init_sign = copysign(1, arl0 - target_arl)

            diff = np.abs(arl0 - target_arl)
            print(diff)

            diff_prop = diff / target_arl

            if method != "ssDMEWMC" and method != "ssMEWMC":
                if diff_prop > 0.50:
                    dp = 0.75
                elif diff_prop > 0.25:
                    dp = 0.5
                elif diff_prop > 0.125:
                    dp = 0.25
                elif diff_prop > 0.0625:
                    dp = 0.125
                elif diff_prop > 0.01:
                    dp = 0.0625
                elif diff_prop > 0.001:
                    dp = 0.01
                else:
                    dp = 0.005
            else:
                if diff_prop > 0.50:
                    dp = 0.75 * h_init
                elif diff_prop > 0.25:
                    dp = 0.5 * h_init
                elif diff_prop > 0.125:
                    dp = 0.25 * h_init
                elif diff_prop > 0.0625:
                    dp = 0.125 * h_init
                elif diff_prop > 0.01:
                    dp = 0.0625 * h_init
                elif diff_prop > 0.001:
                    dp = 0.01 * h_init
                else:
                    dp = 0.005 * h_init

        if c > 0:
            sign = copysign(1, arl0 - target_arl)

            diff = np.abs(arl0 - target_arl)
            print(diff)

            if sign != init_sign:
                dp = dp / 2
                init_sign = sign

        if dp < 0.0000005:
            dp = 0.0000005

        if diff > delta and arl0 > target_arl:
            h_init = h_init - dp
            if h_init in h_list:
                h_init = h_init + (dp / 2)
                dp = dp / 2
            if h_init < 0.001:
                h_init = 0.001
        elif diff > delta and arl0 < target_arl:
            h_init = h_init + dp
            if h_init in h_list:
                h_init = h_init - (dp / 2)
                dp = dp / 2
            if h_init < 0.001:
                h_init = 0.001
        else:
            h_final = h_init
            break

        # an oscillation problem that can't be worked around at this sample size
        if c > 25 and diff < 2 * delta:
            h_final = h_init
            break

        # emergency break clause
        if c > 100:
            h_final = h_init
            break

        h_list.append(h_init)

        c += 1
    return h_final, arl0


# Function to find the out of control run length given an h value, a dimension, a delta, and wait_time
def find_out_of_control_run_length(method, iterations, dimension, smoothing, h_val, delta, change_after):
    print("""Finding the out of control run length for a %s of dimension %d with smoothing parameter %f with a 
          magnitude %f change at time t = %d.""" % (
    method, int(dimension), float(smoothing), float(delta), int(change_after)))

    arl = simulate(method=method, dimension=int(dimension), smoothing=float(smoothing), h=h_val,
                   iterations=int(iterations),
                   out_of_control=True, delta=float(delta), change_after=int(change_after))

    # so the out of control run length is the length of time AFTER the change until it is  detected.
    # The simulate function returns the run length of the whole run. Subtract out the in-control length from the mean.
    arl = arl - (int(dimension) + int(change_after))

    # remove the false positives
    arl = arl[arl > 0]
    arl1 = np.nanmean(arl)

    # print it to the repl
    print(arl1)

    # and return it!
    return arl1

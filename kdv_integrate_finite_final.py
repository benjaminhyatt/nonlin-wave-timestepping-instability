import numpy as np
from sympy.parsing.mathematica import parse_mathematica
from sympy import var
import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger(__name__)
import scipy.integrate as sp

def integrate_finite(alp, dt, el, cinit, LHS_str, RHS_str, scheme, t_stop):
    global alpha 
    global t_step 
    global L 
    global c0 
    global T
    alpha = alp
    t_step = dt
    L = el
    c0 = cinit
    # for rk443
    df = 0.8    

    FWHM = 4*alpha**(1/2)*c0**(-1/2)*np.log(np.sqrt(2)+1)
    T = FWHM/c0
    eps = t_step/T
        
    global sqrt
    global sech
    global tanh
    global sinh
    global cosh
    sqrt = lambda z: np.sqrt(z, dtype = np.float128)
    sech = lambda z: np.cosh(z, dtype = np.float128)**(-1)
    tanh = lambda z: np.tanh(z, dtype = np.float128)
    sinh = lambda z: np.sinh(z, dtype = np.float128)
    cosh = lambda z: np.cosh(z, dtype = np.float128)
    
    def func(t, c):
        numerator = eval(scheme_RHS_expr)
        denominator = eval(scheme_LHS_expr)
        return numerator/denominator
    
    global scheme_LHS_expr 
    global scheme_RHS_expr
    scheme_LHS_expr = LHS_str
    scheme_RHS_expr = RHS_str
    
    scheme_LHS_expr = scheme_LHS_expr.replace('\[Alpha]', 'Alpha').replace('Derivative[1][c][t]','D[c,t]')
    scheme_LHS_expr = parse_mathematica(scheme_LHS_expr)
    #scheme_LHS_expr = str(scheme_LHS_expr).replace('*D(c, t)*','*').replace('Alpha','alpha').replace('c(t)','c')
    scheme_LHS_expr = str(scheme_LHS_expr).replace('*D(c, t)','').replace('Alpha','alpha').replace('c(t)','c')

    scheme_RHS_expr = scheme_RHS_expr.replace('\[Alpha]', 'Alpha')
    scheme_RHS_expr = parse_mathematica(scheme_RHS_expr)
    scheme_RHS_expr = str(scheme_RHS_expr).replace('Alpha','alpha').replace('c(t)','c')
    
    
    solver = sp.RK45
    t_init = 0
    t = t_init
    c = c0*np.ones(1)
    if scheme == 'sbdf1':
        max_dt = t_stop/1e3
    if scheme == 'sbdf2' or scheme == 'rk222':
        max_dt = t_stop/2e3
    if scheme == 'rk443':
        max_dt = (3/2)*(t_stop/2e4)
    print("T_stop", t_stop)
    
    iteration = 0
    
    tslow_list = []
    c_list = []
    tslow_list.append(t)
    c_list.append(c[0])

    while t < t_stop and np.isfinite(c):
        setup = solver(func, t, c, t_stop, max_step = max_dt, atol=1e-12, rtol=1e-8)
        setup.step()
        t = setup.t
        c = setup.y

        tslow_list.append(t)
        c_list.append(c[0])
        iteration +=1
        if iteration % 100 == 0:
            logger.info("Iter=%i, c=%e, tslow=%e" %(iteration, c[0], t))
        # stopping condition for finite blowup schemes
        if scheme != 'rk443' and (np.isclose(c_list[-1], c_list[-2], rtol = 1e-8) or not np.isfinite(c)):
            if scheme != 'sbdf1':
                logger.info("Breaking at iter=%i, c=%e, tslow=%e, t=%e" %(iteration, c[0], t, t/eps**3))
            else:
                logger.info("Breaking at iter=%i, c=%e, tslow=%e, t=%e" %(iteration, c[0], t, t/eps))
            # remove the "repeated" last indices, since it's time will artifically be one timestep late
            tslow_list = tslow_list[:-1]
            c_list = c_list[:-1] 
            break
        if scheme == 'rk443' and c_list[-1] < df*c0:
            logger.info("Stopping at iter=%i, c=%e, tslow=%e, t=%e" %(iteration, c[0], t, t/eps**3))
            break
    tslow_list = np.array(tslow_list)
    c_list = np.array(c_list)
    if scheme == 'sbdf1':
        t_list = tslow_list/eps
    if scheme == 'sbdf2' or scheme == 'rk222' or scheme == 'rk443':
        t_list = tslow_list/eps**3
    
    return t_list, c_list
    

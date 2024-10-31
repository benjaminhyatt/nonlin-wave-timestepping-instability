# note: rk222 rhs uses an approximation for u_inf not used by the other schemes - if we like, we can 
# come back and use this approximation elsewhere
# 10/15/24: I think the above statement was copied in error from solvability_expr.py... 
# so trying to ascertain the situation now... 
# but I think what's going on is that this file has expressions which
# ALL use the u_inf approximation, whereas the other file uses it just for rk222 (out of necessity) 
# and the other schemes use the "full" version 
# I just want to check that that is the case before I explain it that way--I know *this* is the 
# notebook we end up using, but I also recall that the u_inf approximation isn't super valid for rk443
# in the super super long time limit (though it is still numerically sound for some very long times)


# source: kdveulerfinal.nbi
# verified both of these expressions use the u_inf approximation
def get_expr_sbdf1():
    
    sbdf1_LHS_expr = "Times[2.25, Power[L, -1], Power[c[t], -0.5], Plus[Times[Power[L, 2], \
    Power[c[t], 1.5], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 4]], Times[32, Power[c0, 0.5], \[Alpha], Plus[-1, \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]], \
    Times[8, Power[c[t], 0.5], Plus[Times[4, \[Alpha]], Times[Power[c0, \
    0.5], L, Power[\[Alpha], 0.5], Power[Sech[Times[0.25, L, Power[\
    \[Alpha], -0.5], Power[c[t], 0.5]]], 2]], Times[-8, \[Alpha], \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]]], \
    Times[-4, L, Power[\[Alpha], 0.5], c[t], Plus[Times[2, \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[-3, Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]]], Power[Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]], 3]]]], Derivative[1][c][t]]"
    
    sbdf1_RHS_expr = "Times[0.0857143, Power[c, 1.5], Power[L, -3], T, Plus[Times[c, L], \
    Times[12, Plus[Times[-1, Power[c, 0.5]], Power[c0, 0.5]], Power[\
    \[Alpha], 0.5]]], Power[\[Alpha], -0.5], Plus[Times[4, c, L, \
    Plus[Times[17, c, L], Times[84, Plus[Times[-1, Power[c, 0.5]], \
    Power[c0, 0.5]], Power[\[Alpha], 0.5]]]], Times[2, Plus[Times[17, \
    Power[c, 2], Power[L, 2]], Times[126, c, Plus[Power[c, 0.5], \
    Times[-1, Power[c0, 0.5]]], L, Power[\[Alpha], 0.5]], Times[2520, \
    Plus[c, Times[-2, Power[c, 0.5], Power[c0, 0.5]], c0], \[Alpha]]], \
    Power[Sech[Times[0.25, Power[c, 0.5], L, Power[\[Alpha], -0.5]]], \
    2]], Times[-9, c, L, Plus[Times[3, c, L], Times[364, Plus[Power[c, \
    0.5], Times[-1, Power[c0, 0.5]]], Power[\[Alpha], 0.5]]], \
    Power[Sech[Times[0.25, Power[c, 0.5], L, Power[\[Alpha], -0.5]]], \
    4]], Times[450, Power[c, 2], Power[L, 2], Power[Sech[Times[0.25, \
    Power[c, 0.5], L, Power[\[Alpha], -0.5]]], 6]]], Tanh[Times[0.25, \
    Power[c, 0.5], L, Power[\[Alpha], -0.5]]]]"
    
    return sbdf1_LHS_expr, sbdf1_RHS_expr


# source: kdvdecsbdf2final.nb
# verified that these are based on the u_inf approximation
def get_expr_sbdf2():
    
    sbdf2_LHS_expr = "Times[2.25, Power[L, -1], Power[c[t], -0.5], Plus[Times[Power[L, 2], \
    Power[c[t], 1.5], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 4]], Times[32, Power[c0, 0.5], \[Alpha], Plus[-1, \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]], \
    Times[8, Power[c[t], 0.5], Plus[Times[4, \[Alpha]], Times[Power[c0, \
    0.5], L, Power[\[Alpha], 0.5], Power[Sech[Times[0.25, L, Power[\
    \[Alpha], -0.5], Power[c[t], 0.5]]], 2]], Times[-8, \[Alpha], \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]]], \
    Times[-4, L, Power[\[Alpha], 0.5], c[t], Plus[Times[2, \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[-3, Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]]], Power[Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]], 3]]]], Derivative[1][c][t]]"

    sbdf2_RHS_expr = "Times[0.0428571, Power[L, -5], Power[T, 3], Power[\[Alpha], -1.5], \
    Power[c[t], 2.5], Power[Plus[Times[12, Power[c0, 0.5], \
    Power[\[Alpha], 0.5]], Times[-12, Power[\[Alpha], 0.5], Power[c[t], \
    0.5]], Times[L, c[t]]], 3], Plus[Times[720, Power[c0, 0.5], L, Power[\
    \[Alpha], 0.5], c[t]], Times[-720, L, Power[\[Alpha], 0.5], \
    Power[c[t], 1.5]], Times[172, Power[L, 2], Power[c[t], 2]], \
    Times[-15120, c0, \[Alpha], Power[Sech[Times[0.25, L, Power[\[Alpha], \
    -0.5], Power[c[t], 0.5]]], 2]], Times[30240, Power[c0, 0.5], \
    \[Alpha], Power[c[t], 0.5], Power[Sech[Times[0.25, L, Power[\[Alpha], \
    -0.5], Power[c[t], 0.5]]], 2]], Times[780, Power[c0, 0.5], L, Power[\
    \[Alpha], 0.5], c[t], Power[Sech[Times[0.25, L, Power[\[Alpha], \
    -0.5], Power[c[t], 0.5]]], 2]], Times[-15120, \[Alpha], c[t], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[-780, L, Power[\[Alpha], 0.5], Power[c[t], 1.5], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[86, Power[L, 2], Power[c[t], 2], Power[Sech[Times[0.25, L, \
    Power[\[Alpha], -0.5], Power[c[t], 0.5]]], 2]], Times[45360, c0, \
    \[Alpha], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 4]], Times[-90720, Power[c0, 0.5], \[Alpha], \
    Power[c[t], 0.5], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 4]], Times[-23040, Power[c0, 0.5], L, Power[\
    \[Alpha], 0.5], c[t], Power[Sech[Times[0.25, L, Power[\[Alpha], \
    -0.5], Power[c[t], 0.5]]], 4]], Times[45360, \[Alpha], c[t], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    4]], Times[23040, L, Power[\[Alpha], 0.5], Power[c[t], 1.5], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    4]], Times[117, Power[L, 2], Power[c[t], 2], Power[Sech[Times[0.25, \
    L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], 4]], Times[45900, \
    Power[c0, 0.5], L, Power[\[Alpha], 0.5], c[t], Power[Sech[Times[0.25, \
    L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], 6]], Times[-45900, L, \
    Power[\[Alpha], 0.5], Power[c[t], 1.5], Power[Sech[Times[0.25, L, \
    Power[\[Alpha], -0.5], Power[c[t], 0.5]]], 6]], Times[-3525, Power[L, \
    2], Power[c[t], 2], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 6]], Times[7350, Power[L, 2], Power[c[t], 2], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    8]]], Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]"
    return sbdf2_LHS_expr, sbdf2_RHS_expr


# source: kdvrk222not1uinf.nb
# uses u_inf approximation
def get_expr_rk222(): 

    rk222_LHS_expr = "Times[2.25, Power[L, -1], Power[c[t], -0.5], Plus[Times[Power[L, 2], \
    Power[c[t], 1.5], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 4]], Times[32, Power[c0, 0.5], \[Alpha], Plus[-1, \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]], \
    Times[8, Power[c[t], 0.5], Plus[Times[4, \[Alpha]], Times[Power[c0, \
    0.5], L, Power[\[Alpha], 0.5], Power[Sech[Times[0.25, L, Power[\
    \[Alpha], -0.5], Power[c[t], 0.5]]], 2]], Times[-8, \[Alpha], \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]]], \
    Times[-4, L, Power[\[Alpha], 0.5], c[t], Plus[Times[2, \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[-3, Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]]], Power[Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]], 3]]]], Derivative[1][c][t]]"

    rk222_RHS_expr = "Times[0.0000749251, Power[L, -5], Power[T, 3], Power[\[Alpha], \
    -1.5], Power[c[t], 2.5], Plus[Times[12, Power[c0, 0.5], \
    Power[\[Alpha], 0.5]], Times[-12, Power[\[Alpha], 0.5], Power[c[t], \
    0.5]], Times[L, c[t]]], Plus[Times[-8, L, c[t], Plus[Times[-1235520, \
    Power[c0, 1.5], Power[\[Alpha], 1.5]], Times[3706560, c0, Power[\
    \[Alpha], 1.5], Power[c[t], 0.5]], Times[20592, Plus[Times[-31.402, \
    c0, L, \[Alpha]], Times[-180, Power[c0, 0.5], Power[\[Alpha], 1.5]]], \
    c[t]], Times[-41184, Plus[Times[-31.402, Power[c0, 0.5], L, \
    \[Alpha]], Times[-30, Power[\[Alpha], 1.5]]], Power[c[t], 1.5]], \
    Times[156, Plus[Times[-778.097, Power[c0, 0.5], Power[L, 2], Power[\
    \[Alpha], 0.5]], Times[132, -31.402, L, \[Alpha]]], Power[c[t], 2]], \
    Times[-156, -778.097, Power[L, 2], Power[\[Alpha], 0.5], Power[c[t], \
    2.5]], Times[-7946.08, Power[L, 3], Power[c[t], 3]]]], Times[-4, \
    Plus[Times[51891840, Power[c0, 2], Power[\[Alpha], 2]], \
    Times[-207567360, Power[c0, 1.5], Power[\[Alpha], 2], Power[c[t], \
    0.5]], Times[617760, Plus[Times[2.59798, Power[c0, 1.5], L, Power[\
    \[Alpha], 1.5]], Times[504, c0, Power[\[Alpha], 2]]], c[t]], \
    Times[-1853280, Plus[Times[2.59798, c0, L, Power[\[Alpha], 1.5]], \
    Times[112, Power[c0, 0.5], Power[\[Alpha], 2]]], Power[c[t], 1.5]], \
    Times[-10296, Plus[Times[61.7737, c0, Power[L, 2]], Times[-180, \
    2.59798, Power[c0, 0.5], L, Power[\[Alpha], 0.5]], Times[-5040, \
    \[Alpha]]], \[Alpha], Power[c[t], 2]], Times[20592, L, \
    Plus[Times[61.7737, Power[c0, 0.5], L], Times[30, -2.59798, Power[\
    \[Alpha], 0.5]]], \[Alpha], Power[c[t], 2.5]], Times[234, \
    Plus[Times[-522.509, Power[c0, 0.5], Power[L, 3], Power[\[Alpha], \
    0.5]], Times[44, -61.7737, Power[L, 2], \[Alpha]]], Power[c[t], 3]], \
    Times[-234, -522.509, Power[L, 3], Power[\[Alpha], 0.5], Power[c[t], \
    3.5]], Times[-7946.08, Power[L, 4], Power[c[t], 4]]], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[-18, Plus[Times[-34594560, Power[c0, 2], Power[\[Alpha], \
    2]], Times[138378240, Power[c0, 1.5], Power[\[Alpha], 2], Power[c[t], \
    0.5]], Times[-823680, Plus[Times[-20.9045, Power[c0, 1.5], L, Power[\
    \[Alpha], 1.5]], Times[252, c0, Power[\[Alpha], 2]]], c[t]], \
    Times[2471040, Plus[Times[-20.9045, c0, L, Power[\[Alpha], 1.5]], \
    Times[56, Power[c0, 0.5], Power[\[Alpha], 2]]], Power[c[t], 1.5]], \
    Times[6864, Plus[Times[415.514, c0, Power[L, 2]], Times[-360, \
    -20.9045, Power[c0, 0.5], L, Power[\[Alpha], 0.5]], Times[-5040, \
    \[Alpha]]], \[Alpha], Power[c[t], 2]], Times[-13728, L, \
    Plus[Times[415.514, Power[c0, 0.5], L], Times[60, 20.9045, Power[\
    \[Alpha], 0.5]]], \[Alpha], Power[c[t], 2.5]], Times[-52, \
    Plus[Times[-298.091, Power[c0, 0.5], Power[L, 3], Power[\[Alpha], \
    0.5]], Times[132, -415.514, Power[L, 2], \[Alpha]]], Power[c[t], 3]], \
    Times[52, -298.091, Power[L, 3], Power[\[Alpha], 0.5], Power[c[t], \
    3.5]], Times[-1348.9, Power[L, 4], Power[c[t], 4]]], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    4]], Times[30, L, c[t], Plus[Times[-741312, -34.9045, Power[c0, 1.5], \
    Power[\[Alpha], 1.5]], Times[2223936, -34.9045, c0, Power[\[Alpha], \
    1.5], Power[c[t], 0.5]], Times[10296, Plus[Times[272.414, c0, L, \
    \[Alpha]], Times[-216, -34.9045, Power[c0, 0.5], Power[\[Alpha], \
    1.5]]], c[t]], Times[-20592, Plus[Times[272.414, Power[c0, 0.5], L, \
    \[Alpha]], Times[36, 34.9045, Power[\[Alpha], 1.5]]], Power[c[t], \
    1.5]], Times[-78, Plus[Times[12817.0, Power[c0, 0.5], Power[L, 2], \
    Power[\[Alpha], 0.5]], Times[132, -272.414, L, \[Alpha]]], \
    Power[c[t], 2]], Times[78, 12817.0, Power[L, 2], Power[\[Alpha], \
    0.5], Power[c[t], 2.5]], Times[-2883.94, Power[L, 3], Power[c[t], \
    3]]], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]], 6]], Times[-105, Power[L, 2], Power[c[t], 2], \
    Plus[Times[20592, -50.8383, c0, \[Alpha]], Times[-41184, -50.8383, \
    Power[c0, 0.5], \[Alpha], Power[c[t], 0.5]], Times[-312, \
    Plus[Times[3327.77, Power[c0, 0.5], L, Power[\[Alpha], 0.5]], \
    Times[66, 50.8383, \[Alpha]]], c[t]], Times[312, 3327.77, L, Power[\
    \[Alpha], 0.5], Power[c[t], 1.5]], Times[38438.2, Power[L, 2], \
    Power[c[t], 2]]], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 8]], Times[2835, Power[L, 3], Power[c[t], 3], \
    Plus[Times[208, -115.356, Power[c0, 0.5], Power[\[Alpha], 0.5]], \
    Times[208, 115.356, Power[\[Alpha], 0.5], Power[c[t], 0.5]], \
    Times[5966.54, L, c[t]]], Power[Sech[Times[0.25, L, Power[\[Alpha], \
    -0.5], Power[c[t], 0.5]]], 10]], Times[-1715175, 7.90021, Power[L, \
    4], Power[c[t], 4], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 12]]], Tanh[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]]]"
    
    return rk222_LHS_expr, rk222_RHS_expr


# source: kdvrk443final.nb
# uses the u_inf approximation (which is good for the range of c we care about)
def get_expr_rk443(): 
    
    rk443_LHS_expr = "Times[2.25, Power[L, -1], Power[c[t], -0.5], Plus[Times[Power[L, 2], \
    Power[c[t], 1.5], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 4]], Times[32, Power[c0, 0.5], \[Alpha], Plus[-1, \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]], \
    Times[8, Power[c[t], 0.5], Plus[Times[4, \[Alpha]], Times[Power[c0, \
    0.5], L, Power[\[Alpha], 0.5], Power[Sech[Times[0.25, L, Power[\
    \[Alpha], -0.5], Power[c[t], 0.5]]], 2]], Times[-8, \[Alpha], \
    Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]]]]], \
    Times[-4, L, Power[\[Alpha], 0.5], c[t], Plus[Times[2, \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[-3, Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]]], Power[Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]], 3]]]], Derivative[1][c][t]]"
    
    rk443_RHS_expr = "Times[-0.00000624376, Power[L, -5], Power[T, 3], Power[\[Alpha], \
    -1.5], Power[c[t], 2.5], Plus[Times[12, Power[c0, 0.5], \
    Power[\[Alpha], 0.5]], Times[-12, Power[\[Alpha], 0.5], Power[c[t], \
    0.5]], Times[L, c[t]]], Plus[Times[16, L, c[t], Plus[Times[3912480, \
    Power[c0, 1.5], Power[\[Alpha], 1.5]], Times[-11737440, c0, Power[\
    \[Alpha], 1.5], Power[c[t], 0.5]], Times[6864, Plus[Times[377, c0, L, \
    \[Alpha]], Times[1710, Power[c0, 0.5], Power[\[Alpha], 1.5]]], c[t]], \
    Times[-13728, Plus[Times[377, Power[c0, 0.5], L, \[Alpha]], \
    Times[285, Power[\[Alpha], 1.5]]], Power[c[t], 1.5]], Times[26, \
    Plus[Times[27619, Power[c0, 0.5], Power[L, 2], Power[\[Alpha], 0.5]], \
    Times[99528, L, \[Alpha]]], Power[c[t], 2]], Times[-718094, Power[L, \
    2], Power[\[Alpha], 0.5], Power[c[t], 2.5]], Times[77069, Power[L, \
    3], Power[c[t], 3]]]], Times[8, Plus[Times[-164324160, Power[c0, 2], \
    Power[\[Alpha], 2]], Times[657296640, Power[c0, 1.5], Power[\[Alpha], \
    2], Power[c[t], 0.5]], Times[3706560, Plus[Times[3, Power[c0, 1.5], \
    L, Power[\[Alpha], 1.5]], Times[-266, c0, Power[\[Alpha], 2]]], \
    c[t]], Times[-1235520, Plus[Times[27, c0, L, Power[\[Alpha], 1.5]], \
    Times[-532, Power[c0, 0.5], Power[\[Alpha], 2]]], Power[c[t], 1.5]], \
    Times[15444, Plus[Times[117, c0, Power[L, 2]], Times[2160, Power[c0, \
    0.5], L, Power[\[Alpha], 0.5]], Times[-10640, \[Alpha]]], \[Alpha], \
    Power[c[t], 2]], Times[-277992, L, Plus[Times[13, Power[c0, 0.5], L], \
    Times[40, Power[\[Alpha], 0.5]]], \[Alpha], Power[c[t], 2.5]], \
    Times[52, Plus[Times[13232, Power[c0, 0.5], Power[L, 3], Power[\
    \[Alpha], 0.5]], Times[34749, Power[L, 2], \[Alpha]]], Power[c[t], \
    3]], Times[-688064, Power[L, 3], Power[\[Alpha], 0.5], Power[c[t], \
    3.5]], Times[77069, Power[L, 4], Power[c[t], 4]]], \
    Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], 0.5]]], \
    2]], Times[24, Plus[Times[164324160, Power[c0, 2], Power[\[Alpha], \
    2]], Times[-657296640, Power[c0, 1.5], Power[\[Alpha], 2], \
    Power[c[t], 0.5]], Times[-926640, Plus[Times[179, Power[c0, 1.5], L, \
    Power[\[Alpha], 1.5]], Times[-1064, c0, Power[\[Alpha], 2]]], c[t]], \
    Times[308880, Plus[Times[1611, c0, L, Power[\[Alpha], 1.5]], \
    Times[-2128, Power[c0, 0.5], Power[\[Alpha], 2]]], Power[c[t], 1.5]], \
    Times[5148, \[Alpha], Plus[Times[8479, c0, Power[L, 2]], \
    Times[-96660, Power[c0, 0.5], L, Power[\[Alpha], 0.5]], Times[31920, \
    \[Alpha]]], Power[c[t], 2]], Times[-10296, L, Plus[Times[8479, \
    Power[c0, 0.5], L], Times[-16110, Power[\[Alpha], 0.5]]], \[Alpha], \
    Power[c[t], 2.5]], Times[13, Plus[Times[299287, Power[c0, 0.5], \
    Power[L, 3], Power[\[Alpha], 0.5]], Times[3357684, Power[L, 2], \
    \[Alpha]]], Power[c[t], 3]], Times[-3890731, Power[L, 3], Power[\
    \[Alpha], 0.5], Power[c[t], 3.5]], Times[18016, Power[L, 4], \
    Power[c[t], 4]]], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 4]], Times[30, L, c[t], Plus[Times[264771936, \
    Power[c0, 1.5], Power[\[Alpha], 1.5]], Times[-794315808, c0, Power[\
    \[Alpha], 1.5], Power[c[t], 0.5]], Times[-61776, Plus[Times[4189, c0, \
    L, \[Alpha]], Times[-12858, Power[c0, 0.5], Power[\[Alpha], 1.5]]], \
    c[t]], Times[123552, Plus[Times[4189, Power[c0, 0.5], L, \[Alpha]], \
    Times[-2143, Power[\[Alpha], 1.5]]], Power[c[t], 1.5]], Times[-78, \
    Plus[Times[270913, Power[c0, 0.5], Power[L, 2], Power[\[Alpha], \
    0.5]], Times[3317688, L, \[Alpha]]], Power[c[t], 2]], Times[21131214, \
    Power[L, 2], Power[\[Alpha], 0.5], Power[c[t], 2.5]], Times[512177, \
    Power[L, 3], Power[c[t], 3]]], Power[Sech[Times[0.25, L, Power[\
    \[Alpha], -0.5], Power[c[t], 0.5]]], 6]], Times[-105, Power[L, 2], \
    Power[c[t], 2], Plus[Times[-92849328, c0, \[Alpha]], Times[185698656, \
    Power[c0, 0.5], \[Alpha], Power[c[t], 0.5]], Times[-312, \
    Plus[Times[1715, Power[c0, 0.5], L, Power[\[Alpha], 0.5]], \
    Times[297594, \[Alpha]]], c[t]], Times[535080, L, Power[\[Alpha], \
    0.5], Power[c[t], 1.5]], Times[1449353, Power[L, 2], Power[c[t], \
    2]]], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]], 8]], Times[1890, Power[L, 3], Power[c[t], 3], \
    Plus[Times[622492, Power[c0, 0.5], Power[\[Alpha], 0.5]], \
    Times[-622492, Power[\[Alpha], 0.5], Power[c[t], 0.5]], Times[150691, \
    L, c[t]]], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], \
    Power[c[t], 0.5]]], 10]], Times[-112058100, Power[L, 4], Power[c[t], \
    4], Power[Sech[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]], 12]]], Tanh[Times[0.25, L, Power[\[Alpha], -0.5], Power[c[t], \
    0.5]]]]"
    
    return rk443_LHS_expr, rk443_RHS_expr


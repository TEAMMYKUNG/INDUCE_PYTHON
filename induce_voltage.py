from icecream import ic
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations


def pol2cart(r, theta):  # define function polar form <--> cartesian
    z = r * np.exp(1j * theta)
    x, y = z.real, z.imag
    return complex(x, y)


def cart2pol(x, y):
    z = x + y * 1j
    r, theta = np.abs(z), np.angle(z) * 57.2958
    # 1rad = 57.2958 deg
    return r, theta


def cal_distance(Xa, Xb, Ya, Yb):  # *function find Distance between 2 point
    distance = np.sqrt((Xa - Xb) ** 2 + (Ya - Yb) ** 2)
    return distance


def calculate(xp, yp, obj_size):
    xp = -xp
    ###############################################
    if case == 1:  # case 1: 115kv
        (xa, ya) = (0, 15.1)
        (xb, yb) = (4, 15.1)
        (xc, yc) = (4, 12.1)
        (r_a, theta_a) = (115, 0)
        (r_b, theta_b) = (115, -120)
        (r_c, theta_c) = (115, 120)
        conductor = ('a', 'b', 'c')
        comb = combinations([conductor[i] for i in range(len(conductor))], 2)
        distant, distantp, distantD, vphase, iphase = {}, {}, {}, {}, {}

        for i in list(comb):
            key = str(str(i[0]) + str(i[1]))
            distant[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]], vars()['y' + i[1]])
            distantp[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]],
                                         - vars()['y' + i[1]])
        for i in conductor:
            key = str('p' + i)
            distant[key] = cal_distance(xp, vars()['x' + i], yp, vars()['y' + i])
            distantp[key] = cal_distance(xp, vars()['x' + i], yp, - vars()['y' + i])
            distantD[key] = cal_distance(xp, vars()['x' + i], 0, vars()['y' + i])
            vphase[i] = pol2cart(vars()['r_' + i], np.radians(vars()['theta_' + i]))
            iphase[i] = pol2cart(line_current_115, np.radians(vars()['theta_' + i]))

        v_complex = np.array([[vphase['a'] / np.sqrt(3)], [vphase['b'] / np.sqrt(3)], [vphase['c'] / np.sqrt(3)]])
        Matrix = np.array(
            [[np.log(2 * ya / r_115), np.log(distantp['ab'] / distant['ab']), np.log(distantp['ac'] / distant['ac'])],
             [np.log(distantp['ab'] / distant['ab']), np.log(2 * yb / r_115), np.log(distantp['bc'] / distant['bc'])],
             [np.log(distantp['ac'] / distant['ac']), np.log(distantp['bc'] / distant['bc']), np.log(2 * yc / r_115)]])
        q_cart = 2 * np.pi * EPSILON_0 * np.matmul(np.linalg.inv(Matrix), v_complex)
        vpe = ((q_cart[0] * np.log(distantp['pa'] / distant['pa'])) + (
                q_cart[1] * np.log(distantp['pb'] / distant['pb'])) + (
                       q_cart[2] * np.log(distantp['pc'] / distant['pc']))) / (2 * np.pi * EPSILON_0)
        i_complex = np.array([[iphase['a']], [iphase['b']], [iphase['c']]])
        SuperPosition = np.array([[np.log(distantD['pa'] / distant['pa']), np.log(distantD['pb'] / distant['pb']),
                                   np.log(distantD['pc'] / distant['pc'])]])
        matrix2 = np.matmul(SuperPosition, i_complex)
        ep = 2 * (10 ** -7) * 100 * np.pi * matrix2
        vpm = ep * obj_size
        vp = vpm + vpe
        (VP, VI) = cart2pol(np.real(vp), np.imag(vp))
        return round(VP[0][0], 2)

    ############################################################################################
    if case == 2:  # case 2: 115kv+22kv
        (xa, ya) = (4, 18.8)
        (xb, yb) = (4, 16.3)
        (xc, yc) = (4, 13.8)
        (xd, yd) = (0.85, 13.28)
        (xe, ye) = (2.6, 13.28)
        (xf, yf) = (3.15, 13.28)
        (r_a, theta_a) = (115, 0)
        (r_b, theta_b) = (115, -120)
        (r_c, theta_c) = (115, 120)
        (r_d, theta_d) = (22, 0)
        (r_e, theta_e) = (22, -120)
        (r_f, theta_f) = (22, 120)
        conductor = ('a', 'b', 'c', 'd', 'e', 'f')
        comb = combinations([conductor[i] for i in range(len(conductor))], 2)
        distant, distantp, distantD, vphase, iphase = {}, {}, {}, {}, {}
        for i in list(comb):
            key = str(str(i[0]) + str(i[1]))
            distant[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]], vars()['y' + i[1]])
            distantp[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]],
                                         - vars()['y' + i[1]])
        for i in conductor:
            key = str('p' + i)
            distant[key] = cal_distance(xp, vars()['x' + i], yp, vars()['y' + i])
            distantp[key] = cal_distance(xp, vars()['x' + i], yp, - vars()['y' + i])
            distantD[key] = cal_distance(xp, vars()['x' + i], 0, vars()['y' + i])
            vphase[i] = pol2cart(vars()['r_' + i], np.radians(vars()['theta_' + i]))
            if i in ('a', 'b', 'c'):
                iphase[i] = pol2cart(line_current_115, np.radians(vars()['theta_' + i]))
            if i in ('d', 'e', 'f'):
                iphase[i] = pol2cart(line_current_22, np.radians(vars()['theta_' + i]))

        v_complex = np.array([[vphase['a'] / np.sqrt(3)], [vphase['b'] / np.sqrt(3)], [vphase['c'] / np.sqrt(3)],
                              [vphase['d'] / np.sqrt(3)], [vphase['e'] / np.sqrt(3)], [vphase['f'] / np.sqrt(3)]])
        Matrix = np.array([[np.log(2 * ya / r_115), np.log(distantp['ab'] / distant['ab']),
                            np.log(distantp['ac'] / distant['ac']), np.log(distantp['ad'] / distant['ad']),
                            np.log(distantp['ae'] / distant['ae']), np.log(distantp['af'] / distant['af'])],
                           [np.log(distantp['ab'] / distant['ab']), np.log(2 * yb / r_115),
                            np.log(distantp['bc'] / distant['bc']), np.log(distantp['bd'] / distant['bd']),
                            np.log(distantp['be'] / distant['be']), np.log(distantp['bf'] / distant['bf'])],
                           [np.log(distantp['ac'] / distant['ac']), np.log(distantp['bc'] / distant['bc']),
                            np.log(2 * yc / r_115), np.log(distantp['cd'] / distant['cd']),
                            np.log(distantp['ce'] / distant['ce']), np.log(distantp['cf'] / distant['cf'])],
                           [np.log(distantp['ad'] / distant['ad']), np.log(distantp['bd'] / distant['bd']),
                            np.log(distantp['cd'] / distant['cd']), np.log(2 * yd / r_22),
                            np.log(distantp['de'] / distant['de']), np.log(distantp['df'] / distant['df'])],
                           [np.log(distantp['ae'] / distant['ae']), np.log(distantp['be'] / distant['be']),
                            np.log(distantp['ce'] / distant['ce']), np.log(distantp['de'] / distant['de']),
                            np.log(2 * ye / r_22), np.log(distantp['ef'] / distant['ef'])],
                           [np.log(distantp['af'] / distant['af']), np.log(distantp['bf'] / distant['bf']),
                            np.log(distantp['cf'] / distant['cf']), np.log(distantp['df'] / distant['df']),
                            np.log(distantp['ef'] / distant['ef']), np.log(2 * yf / r_22)]])
        q_cart = 2 * np.pi * EPSILON_0 * np.matmul(np.linalg.inv(Matrix), v_complex)
        vpe = ((q_cart[0] * np.log(distantp['pa'] / distant['pa'])) + (
                q_cart[1] * np.log(distantp['pb'] / distant['pb'])) + (
                       q_cart[2] * np.log(distantp['pc'] / distant['pc'])) + (
                       q_cart[3] * np.log(distantp['pd'] / distant['pd'])) + (
                       q_cart[4] * np.log(distantp['pe'] / distant['pe'])) + (
                       q_cart[5] * np.log(distantp['pf'] / distant['pf']))) / (2 * np.pi * EPSILON_0)
        i_complex = np.array([[iphase['a']], [iphase['b']], [iphase['c']], [iphase['d']], [iphase['e']], [iphase['f']]])
        SuperPosition = np.array([[np.log(distantD['pa'] / distant['pa']), np.log(distantD['pb'] / distant['pb']),
                                   np.log(distantD['pc'] / distant['pc']), np.log(distantD['pd'] / distant['pd']),
                                   np.log(distantD['pe'] / distant['pe']), np.log(distantD['pf'] / distant['pf'])]])
        matrix2 = np.matmul(SuperPosition, i_complex)
        ep = 2 * (10 ** -7) * 100 * np.pi * matrix2
        vpm = ep * obj_size
        vp = vpm + vpe
        (VP, VI) = cart2pol(np.real(vp), np.imag(vp))
        return round(VP[0][0], 2)

    ############################################################################################
    if case == 3:  # case 3: 115kv Double Circuit + 22kv
        (xa, ya) = (4, 18.8)
        (xb, yb) = (4, 16.3)
        (xc, yc) = (4, 13.8)
        (xd, yd) = (0.85, 13.28)
        (xe, ye) = (2.6, 13.28)
        (xf, yf) = (3.15, 13.28)
        (xg, yg) = (0, 18.8)
        (xh, yh) = (0, 17.3)
        (xi, yi) = (0, 13.8)
        (r_a, theta_a) = (115, 0)
        (r_b, theta_b) = (115, -120)
        (r_c, theta_c) = (115, 120)
        (r_d, theta_d) = (22, 0)
        (r_e, theta_e) = (22, -120)
        (r_f, theta_f) = (22, 120)
        (r_g, theta_g) = (115, 0)
        (r_h, theta_h) = (115, -120)
        (r_i, theta_i) = (115, 120)
        conductor = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')
        comb = combinations([conductor[i] for i in range(len(conductor))], 2)
        distant, distantp, distantD, vphase, iphase = {}, {}, {}, {}, {}
        for i in list(comb):
            key = str(str(i[0]) + str(i[1]))
            distant[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]], vars()['y' + i[1]])
            distantp[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]], - vars()['y' + i[1]])
        for i in conductor:
            key = str('p' + i)
            distant[key] = cal_distance(xp, vars()['x' + i], yp, vars()['y' + i])
            distantp[key] = cal_distance(xp, vars()['x' + i], yp, - vars()['y' + i])
            distantD[key] = cal_distance(xp, vars()['x' + i], 0, vars()['y' + i])
            vphase[i] = pol2cart(vars()['r_' + i], np.radians(vars()['theta_' + i]))
            if i in ('a', 'b', 'c', 'g', 'h', 'i'):
                iphase[i] = pol2cart(line_current_115, np.radians(vars()['theta_' + i]))
            if i in ('d', 'e', 'f'):
                iphase[i] = pol2cart(line_current_22, np.radians(vars()['theta_' + i]))

        v_complex = np.array([[vphase['a'] / np.sqrt(3)], [vphase['b'] / np.sqrt(3)], [vphase['c'] / np.sqrt(3)],
                              [vphase['d'] / np.sqrt(3)], [vphase['e'] / np.sqrt(3)], [vphase['f'] / np.sqrt(3)],
                              [vphase['g'] / np.sqrt(3)], [vphase['h'] / np.sqrt(3)], [vphase['i'] / np.sqrt(3)]])
        Matrix = np.array([[np.log(2 * ya / r_115), np.log(distantp['ab'] / distant['ab']),
                            np.log(distantp['ac'] / distant['ac']), np.log(distantp['ad'] / distant['ad']),
                            np.log(distantp['ae'] / distant['ae']), np.log(distantp['af'] / distant['af']),
                            np.log(distantp['ag'] / distant['ag']), np.log(distantp['ah'] / distant['ah']),
                            np.log(distantp['ai'] / distant['ai'])],
                            [np.log(distantp['ab'] / distant['ab']), np.log(2 * yb / r_115),
                            np.log(distantp['bc'] / distant['bc']), np.log(distantp['bd'] / distant['bd']),
                            np.log(distantp['be'] / distant['be']), np.log(distantp['bf'] / distant['bf']),
                            np.log(distantp['bg'] / distant['bg']), np.log(distantp['bh'] / distant['bh']),
                            np.log(distantp['bi'] / distant['bi'])],
                            [np.log(distantp['ac'] / distant['ac']), np.log(distantp['bc'] / distant['bc']),
                            np.log(2 * yc / r_115), np.log(distantp['cd'] / distant['cd']),
                            np.log(distantp['ce'] / distant['ce']), np.log(distantp['cf'] / distant['cf']),
                            np.log(distantp['cg'] / distant['cg']), np.log(distantp['ch'] / distant['ch']),
                            np.log(distantp['ci'] / distant['ci'])],
                            [np.log(distantp['ad'] / distant['ad']), np.log(distantp['bd'] / distant['bd']),
                            np.log(distantp['cd'] / distant['cd']), np.log(2 * yd / r_22),
                            np.log(distantp['de'] / distant['de']), np.log(distantp['df'] / distant['df']),
                            np.log(distantp['dg'] / distant['dg']), np.log(distantp['dh'] / distant['dh']),
                            np.log(distantp['di'] / distant['di'])],
                            [np.log(distantp['ae'] / distant['ae']), np.log(distantp['be'] / distant['be']),
                            np.log(distantp['ce'] / distant['ce']), np.log(distantp['de'] / distant['de']),
                            np.log(2 * ye / r_22), np.log(distantp['ef'] / distant['ef']),
                            np.log(distantp['eg'] / distant['eg']), np.log(distantp['eh'] / distant['eh']),
                            np.log(distantp['ei'] / distant['ei'])],
                            [np.log(distantp['af'] / distant['af']), np.log(distantp['bf'] / distant['bf']),
                            np.log(distantp['cf'] / distant['cf']), np.log(distantp['df'] / distant['df']),
                            np.log(distantp['ef'] / distant['ef']), np.log(2 * yf / r_22),
                            np.log(distantp['fg'] / distant['fg']), np.log(distantp['fh'] / distant['fh']),
                            np.log(distantp['fi'] / distant['fi'])],
                            [np.log(distantp['ag'] / distant['ag']), np.log(distantp['bg'] / distant['bg']),
                            np.log(distantp['cg'] / distant['cg']), np.log(distantp['dg'] / distant['dg']),
                            np.log(distantp['eg'] / distant['eg']), np.log(distantp['fg'] / distant['fg']),
                            np.log(2 * yg / r_115), np.log(distantp['gh'] / distant['gh']),
                            np.log(distantp['gi'] / distant['gi'])],
                            [np.log(distantp['ah'] / distant['ah']), np.log(distantp['bh'] / distant['bh']),
                            np.log(distantp['ch'] / distant['ch']), np.log(distantp['dh'] / distant['dh']),
                            np.log(distantp['eh'] / distant['eh']), np.log(distantp['fh'] / distant['fh']),
                            np.log(distantp['gh'] / distant['gh']),np.log(2 * yh / r_115),
                            np.log(distantp['hi'] / distant['hi'])],
                            [np.log(distantp['ai'] / distant['ai']), np.log(distantp['bi'] / distant['bi']),
                            np.log(distantp['ci'] / distant['ci']), np.log(distantp['di'] / distant['di']),
                            np.log(distantp['ei'] / distant['ei']), np.log(distantp['fi'] / distant['fi']),
                            np.log(distantp['gi'] / distant['gi']),np.log(distantp['hi'] / distant['hi']),
                            np.log(2 * yi / r_115)]])

        q_cart = 2 * np.pi * EPSILON_0 * np.matmul(np.linalg.inv(Matrix), v_complex)
        vpe = ((q_cart[0] * np.log(distantp['pa'] / distant['pa'])) +
               (q_cart[1] * np.log(distantp['pb'] / distant['pb'])) +
               (q_cart[2] * np.log(distantp['pc'] / distant['pc'])) +
               (q_cart[3] * np.log(distantp['pd'] / distant['pd'])) +
               (q_cart[4] * np.log(distantp['pe'] / distant['pe'])) +
               (q_cart[5] * np.log(distantp['pf'] / distant['pf'])) +
               (q_cart[5] * np.log(distantp['pg'] / distant['pg'])) +
               (q_cart[5] * np.log(distantp['ph'] / distant['ph'])) +
               (q_cart[5] * np.log(distantp['pi'] / distant['pi']))) / (2 * np.pi * EPSILON_0)
        i_complex = np.array([[iphase['a']], [iphase['b']], [iphase['c']], [iphase['d']], [iphase['e']], [iphase['f']], [iphase['g']], [iphase['h']], [iphase['i']]])
        SuperPosition = np.array([[np.log(distantD['pa'] / distant['pa']), np.log(distantD['pb'] / distant['pb']),
                                   np.log(distantD['pc'] / distant['pc']), np.log(distantD['pd'] / distant['pd']),
                                   np.log(distantD['pe'] / distant['pe']), np.log(distantD['pf'] / distant['pf']),
                                   np.log(distantD['pg'] / distant['pg']), np.log(distantD['ph'] / distant['ph']),
                                   np.log(distantD['pi'] / distant['pi'])]])
        matrix2 = np.matmul(SuperPosition, i_complex)
        ep = 2 * (10 ** -7) * 100 * np.pi * matrix2
        vpm = ep * obj_size
        vp = vpm + vpe
        (VP, VI) = cart2pol(np.real(vp), np.imag(vp))
        return round(VP[0][0], 2)

    ############################################################################################
    if case == 4:  # case 4: 115kv Bundle
        (xa, ya) = (1.95, 18.8)
        (xb, yb) = (2.15, 18.8)
        (xc, yc) = (1.95, 16.3)
        (xd, yd) = (2.15, 16.3)
        (xe, ye) = (1.95, 13.8)
        (xf, yf) = (2.15, 13.8)
        (r_a, theta_a) = (115, 0)
        (r_b, theta_b) = (115, 0)
        (r_c, theta_c) = (115, 120)
        (r_d, theta_d) = (115, 120)
        (r_e, theta_e) = (115, -120)
        (r_f, theta_f) = (115, -120)
        conductor = ('a', 'b', 'c', 'd', 'e', 'f')
        comb = combinations([conductor[i] for i in range(len(conductor))], 2)
        distant, distantp, distantD, vphase, iphase = {}, {}, {}, {}, {}
        for i in list(comb):
            key = str(str(i[0]) + str(i[1]))
            distant[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]], vars()['y' + i[1]])
            distantp[key] = cal_distance(vars()['x' + i[0]], vars()['x' + i[1]], vars()['y' + i[0]],
                                         - vars()['y' + i[1]])
        for i in conductor:
            key = str('p' + i)
            distant[key] = cal_distance(xp, vars()['x' + i], yp, vars()['y' + i])
            distantp[key] = cal_distance(xp, vars()['x' + i], yp, - vars()['y' + i])
            distantD[key] = cal_distance(xp, vars()['x' + i], 0, vars()['y' + i])
            vphase[i] = pol2cart(vars()['r_' + i], np.radians(vars()['theta_' + i]))
            iphase[i] = pol2cart(line_current_115, np.radians(vars()['theta_' + i]))

        v_complex = np.array([[vphase['a'] / np.sqrt(3)], [vphase['b'] / np.sqrt(3)], [vphase['c'] / np.sqrt(3)],
                              [vphase['d'] / np.sqrt(3)], [vphase['e'] / np.sqrt(3)], [vphase['f'] / np.sqrt(3)]])
        Matrix = np.array([[np.log(2 * ya / gmr_115), np.log(distantp['ab'] / distant['ab']),
                            np.log(distantp['ac'] / distant['ac']), np.log(distantp['ad'] / distant['ad']),
                            np.log(distantp['ae'] / distant['ae']), np.log(distantp['af'] / distant['af'])],
                           [np.log(distantp['ab'] / distant['ab']), np.log(2 * yb / gmr_115),
                            np.log(distantp['bc'] / distant['bc']), np.log(distantp['bd'] / distant['bd']),
                            np.log(distantp['be'] / distant['be']), np.log(distantp['bf'] / distant['bf'])],
                           [np.log(distantp['ac'] / distant['ac']), np.log(distantp['bc'] / distant['bc']),
                            np.log(2 * yc / gmr_115), np.log(distantp['cd'] / distant['cd']),
                            np.log(distantp['ce'] / distant['ce']), np.log(distantp['cf'] / distant['cf'])],
                           [np.log(distantp['ad'] / distant['ad']), np.log(distantp['bd'] / distant['bd']),
                            np.log(distantp['cd'] / distant['cd']), np.log(2 * yd / gmr_115),
                            np.log(distantp['de'] / distant['de']), np.log(distantp['df'] / distant['df'])],
                           [np.log(distantp['ae'] / distant['ae']), np.log(distantp['be'] / distant['be']),
                            np.log(distantp['ce'] / distant['ce']), np.log(distantp['de'] / distant['de']),
                            np.log(2 * ye / gmr_115), np.log(distantp['ef'] / distant['ef'])],
                           [np.log(distantp['af'] / distant['af']), np.log(distantp['bf'] / distant['bf']),
                            np.log(distantp['cf'] / distant['cf']), np.log(distantp['df'] / distant['df']),
                            np.log(distantp['ef'] / distant['ef']), np.log(2 * yf / gmr_115)]])
        q_cart = 2 * np.pi * EPSILON_0 * np.matmul(np.linalg.inv(Matrix), v_complex)
        vpe = ((q_cart[0] * np.log(distantp['pa'] / distant['pa'])) + (
                q_cart[1] * np.log(distantp['pb'] / distant['pb'])) + (
                q_cart[2] * np.log(distantp['pc'] / distant['pc'])) + (
                q_cart[3] * np.log(distantp['pd'] / distant['pd'])) + (
                q_cart[4] * np.log(distantp['pe'] / distant['pe'])) + (
                q_cart[5] * np.log(distantp['pf'] / distant['pf']))) / (2 * np.pi * EPSILON_0)
        i_complex = np.array([[iphase['a']], [iphase['b']], [iphase['c']], [iphase['d']], [iphase['e']], [iphase['f']]])
        SuperPosition = np.array([[np.log(distantD['pa'] / distant['pa']), np.log(distantD['pb'] / distant['pb']),
                                   np.log(distantD['pc'] / distant['pc']), np.log(distantD['pd'] / distant['pd']),
                                   np.log(distantD['pe'] / distant['pe']), np.log(distantD['pf'] / distant['pf'])]])
        matrix2 = np.matmul(SuperPosition, i_complex)
        ep = 2 * (10 ** -7) * 100 * np.pi * matrix2
        vpm = ep * obj_size
        vp = vpm + vpe
        (VP, VI) = cart2pol(np.real(vp), np.imag(vp))
        return round(VP[0][0], 2)

    else:
        print("Don't have that case")
        exit()


# define Gobal var
r_115 = 12.8 #0.012825
r_22 = 7.9 #0.00799
gmr_115 = 50.69 #0.05069
line_current_115 = 1606.539
line_current_22 = 1606.539  #############
EPSILON_0 = 8.854 * pow(10, -12)
v_safe = 50
con = 'y'

# Start
while con == 'y' or con == 'Y':
    case = int(input('Enter Case: '))
    obj_size = int(input('Enter Size: '))
    inter_h = float(input('Enter interested height : '))

    # Graph on interested High
    m_xp = []
    v_induce = []
    for Distance_x in np.linspace(0.1, 10, 100):
        Distance_x = round(Distance_x, 1)
        m_xp.append(Distance_x)
        v_induce.append(calculate(Distance_x, inter_h, obj_size))
    data = {'distance': m_xp, 'Induced Voltage': v_induce}
    df = pd.DataFrame.from_dict(data)
    df.set_index('distance', inplace=True)
    df.plot(figsize=(10, 8), ylabel='Induce Voltage(V)', xlabel='Distance(m)', title='Induce Voltage & Distance',
            grid=True, xlim=0, ylim=0)

    # heatmap
    x_ax = np.arange(0.1, 5, 0.1)
    y_ax = np.arange(0.1, 20, 0.1)
    x_ax = np.round(x_ax, 2)
    y_ax = np.round(y_ax, 2)
    z_ax = np.zeros(shape=(len(y_ax), len(x_ax)))
    z_con = np.zeros(shape=(len(y_ax), len(x_ax)))
    cx = -1
    for x_vax in x_ax:
        cx += 1
        cy = -1
        for y_vax in y_ax:
            cy += 1
            z_ax[cy][cx] = calculate(x_vax, y_vax, obj_size)
            if calculate(x_vax, y_vax, obj_size) >= v_safe:
                z_con[cy][cx] = 1
            else:
                z_con[cy][cx] = 0
    # Normal Heatmap
    df = pd.DataFrame(data=z_ax, columns=x_ax, index=y_ax)
    plt.figure(figsize=(16, 9))
    heatmap = sns.heatmap(df, cbar_kws={'label': 'Induce Voltage(V)'}, cmap='OrRd')
    heatmap.invert_yaxis()
    heatmap.invert_xaxis()
    heatmap.set(xlabel='Distance(m)', ylabel='High(m)')

    # safe zone Heatmap
    df = pd.DataFrame(data=z_con, columns=x_ax, index=y_ax)
    plt.figure(figsize=(16, 9))
    heatmap = sns.heatmap(df, cmap='OrRd', cbar=False)
    heatmap.invert_yaxis()
    heatmap.invert_xaxis()
    heatmap.set(xlabel='Distance(m)', ylabel='High(m)')
    plt.show()

    print('Calculate Again ?')
    con = str(input('Enter y or n: '))
else:
    exit()

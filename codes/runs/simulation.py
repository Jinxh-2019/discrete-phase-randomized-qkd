from codes.runs.startpoint import Xs, Lgts
from codes.runs.sql_exe import create_table, write_records, connect,get_initials,put_initials
from codes.utils.optimize import global_optimize, optimize
from codes.utils.keyrate import keyrate
from codes.utils.sampling import sampling
from codes.utils.search_for_initial_point import search_for_initial_point
from sqlite3 import Cursor
from numpy import Inf


def generate_testitems(kwargs):
    res = []
    for mode in kwargs['modes']:
        if mode in {'DP'}:
            for lemma_ver in kwargs['lemmas']:
                for logNtot in kwargs['logNtots']:
                    table_name = lemma_ver + '_' + mode + '10_%d' % logNtot
                    res.append((mode, table_name, logNtot,lemma_ver))
        elif mode in {'CPFK'}:
            for logNtot in kwargs['logNtots']:
                table_name = mode + '10_%d' % logNtot
                res.append((mode, table_name, logNtot))
        elif mode in{'CP_Ma','DP_Cao'}:
            res.append((mode,mode))
    return res


def set_through_testitem(kwargs, testitem):
    mode = testitem[0]
    kwargs['mode'] = mode
    if mode in {'DP'}:
        lemma_ver = testitem[3]
        table_name = testitem[1]
        logNtot = testitem[2]
        table_name = lemma_ver + '_' + mode +'10_%d' % logNtot
        kwargs['table_name'] = table_name
        kwargs['logNtot'] = logNtot
        kwargs['Ntot'] = 10**logNtot
        kwargs['mode_Ntot'] = mode + '10_%d' % logNtot
        kwargs['lemma'] = lemma_ver
    elif mode in {'CPFK'}:
        table_name = testitem[1]
        logNtot = testitem[2]
        table_name = mode +'10_%d' % logNtot
        kwargs['table_name'] = table_name
        kwargs['logNtot'] = logNtot
        kwargs['Ntot'] = 10**logNtot
        kwargs['mode_Ntot'] = mode + '10_%d' % logNtot
    elif mode in {'DP_Cao','CP_Ma'}:
        table_name = testitem[1]
        kwargs['table_name'] = table_name
        kwargs['mode_Ntot'] = mode

        # kwargs['logNtot'] = None
        # kwargs['Ntot'] = None
    return table_name

def str_to_float(X):
    items = X[1:-2].split(',')
    ans = []
    for item in items:
        ans +=[float(item)]
    return ans
def simulator(if_optimize: bool, kwargs):
    cur, conn = connect(kwargs['conf'], kwargs)
    testitems = generate_testitems(kwargs)
    initials = get_initials()
    initials = dict(initials)
    
    for testitem in testitems:
        table_name = set_through_testitem(kwargs, testitem)
        if table_name in initials:
            X = str_to_float(initials[table_name])
        else:
            X = search_for_initial_point(kwargs)
            put_initials(table_name,X)
        create_table(cur, table_name, kwargs)
        simulator_for_per_testitem(
            cur, table_name, if_optimize, kwargs, X)
        conn.commit()
    conn.close()

def simulator_for_per_point(kwargs,l,xini,if_optimize):
    mode = kwargs['mode']
    if if_optimize:
        if mode in {'DP','CPFK'}:
            x, b = global_optimize(kwargs, l=l, x=xini)
        elif mode in {'DP_Cao', 'CP_Ma'}:
            x, b = optimize(kwargs, l=l, x=xini)
        if x == []:
            print('optimize error')
            key_rate = -Inf
        else:
            key_rate = keyrate(x, l, kwargs)
    else:
        key_rate = keyrate(xini, l, kwargs)
        x = xini
    return list(x),key_rate
    
def simulator_for_per_testitem(cur: Cursor, table_name, if_optimize, kwargs, start_x=[]):
    if start_x == []: x = sampling()
    else: x = start_x
    count_down = 5
    last_eff_l = 0
    last_eff_x = x 
    for l in range(0,1000,10):
        kwargs['Lgt'] = l
        x,key_rate = simulator_for_per_point(kwargs,l,last_eff_x,if_optimize)
        if x != []:
            last_eff_x = x
        else:
            continue
        if key_rate >= 0:
            last_eff_l = l
            write_records(cur, x, l, key_rate, table_name)
            print(table_name, l, key_rate, '\n', x)
            count_down = 5
        else:
            count_down = count_down-1
            if count_down == 0:
                break
    
    for l in range(last_eff_l-10,last_eff_l+10):
        kwargs['Lgt'] = l
        x,key_rate = simulator_for_per_point(kwargs,l,last_eff_x,if_optimize)
        if x != []:
            last_eff_x = x
        else:
            continue
        if key_rate >= 0:
            # last_eff_l = l
            write_records(cur, x, l, key_rate, table_name)
            print(table_name, l, key_rate, '\n', x)
            count_down = 5
        else:
            count_down = count_down-1
            if count_down == 0:
                break
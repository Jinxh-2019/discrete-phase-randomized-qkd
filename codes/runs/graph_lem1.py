import matplotlib.pyplot as plt
import numpy as np
from sqlite3 import Cursor
from codes.runs.sql_exe import search, connect
from codes.runs.simulation import generate_testitems


def draw_single_line(cur: Cursor, table_name, label, discard):
    data = search(cur, table_name)
    X = []
    Y = []
    for p in data:
        X.append(p[1])
        Y.append(p[2])
    discard_list = []
    for i in range(len(X)):
        if X[i] in discard:
            discard_list.append(i)
    for j in range(len(discard_list)-1, -1, -1):
        del X[discard_list[j]]
        del Y[discard_list[j]]
    plt.semilogy(X, Y, label=label)


discards = {'DP10_12': {}, 'DP10_13': {}, 'DP10_14': {},
            'DP10_20': {}, 'DP_Cao': {}, 'CP_Ma': {}}


def output_graphs(kwargs):
    cur, conn = connect(kwargs["conf"], kwargs)
    testitems = generate_testitems(kwargs)
    for item in testitems:
        table_name = item[2]
        draw_single_line(
            cur, table_name, table_name, discards[kwargs["mode"] + '10_%d' % kwargs["logNtot"]])
    conn.close()
    X = np.array(range(320))
    Y = 0.3*10**(-X*0.2/10)
    plt.semilogy(X, Y, label='PLOB linear bound')

    plt.xlabel('Distance (km)')
    plt.ylabel('Keyrate (bits/pulse)')
    plt.legend()
    #plt.savefig('fig3.jpg',dpi=1200,figsize = (24,32))
    plt.show()
    print('end')

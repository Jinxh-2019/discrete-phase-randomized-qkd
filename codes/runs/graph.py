import matplotlib.pyplot as plt
import numpy as np
from sqlite3 import Cursor
from codes.runs.sql_exe import search, connect
from codes.runs.simulation import generate_testitems


def draw_single_line(cur: Cursor, table_name, label, discard={}):
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
    for j in range(len(discard_list)-1,-1,-1):
        del X[discard_list[j]]
        del Y[discard_list[j]]
    if table_name == 'CP_Ma':
        plt.semilogy(X, Y, label=label,linestyle='--')
    elif table_name == 'DP_Cao':
        plt.semilogy(X, Y, label=label,linestyle='-.')
    else:
        plt.semilogy(X, Y, label=label)


discards = {'lemmaA1_DP10_20': {20,30,40,50,60,70,80,90}
            }

def graph_name(table_name):
    name_dict = {'DP_Cao': 'Discrete Phase Randomize asymptotic',
                    'CP_Ma': 'Continuous Randomize Phase asymptotic',
                    'lemmaA1_DP10_12': 'Discrete Phase Randomize with Ntot=10^12',
                    'lemmaA1_DP10_14': 'Discrete Phase Randomize with Ntot=10^14',
                    'lemmaA1_DP10_20': 'Discrete Phase Randomize with Ntot=10^20',
                    'lemma_shan_DP10_20': 'Discrete Phase Randomize with Ntot=10^20,shan\'s method',
                    'lemma_shan_DP10_14': 'Discrete Phase Randomize with Ntot=10^14,shan\'s method',
                    'lemma_shan_DP10_12': 'Discrete Phase Randomize with Ntot=10^12,shan\'s method',
                    }
    if table_name in name_dict:
        return name_dict[table_name]
    else: return table_name
def output_graphs(kwargs):
    plt.figure(figsize=(10,8))
    cur, conn = connect(kwargs['conf'],kwargs)
    # test_items = generate_testitems(kwargs)
    cur.execute(
        " SELECT name from sqlite_master WHERE type='table'")
    table_names = cur.fetchall()
    for wrapped_table_name in table_names:
        table_name = wrapped_table_name[0]
        if not table_name in {'CP_Ma','DP_Cao','lemmaA1_DP10_12','lemmaA1_DP10_14','lemmaA1_DP10_20'}: continue
        if table_name in discards: discard = discards[table_name]
        else: discard = {}
        draw_single_line(
            cur, table_name, graph_name(table_name), discard)
    conn.close()
    # X = np.array(range(320))
    # Y = 0.3*10**(-X*0.2/10)
    # plt.semilogy(X, Y, label='PLOB linear bound')

    plt.xlabel('Distance (km)')
    plt.ylabel('Keyrate (bits/pulse)')
    plt.legend()
    plt.savefig('fig3.eps',dpi=1200)
    # plt.show()
    print('end')

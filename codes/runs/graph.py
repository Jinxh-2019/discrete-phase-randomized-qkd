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
    plt.semilogy(X, Y, label=label)


discards = {'CP_Ma': {251,252,253,254,255,256,257,258,259,261,262,263,264,265,266,267,268,269}
            }

def graph_name(table_name):
    name_dict = {'DP_Cao': 'Discrete Phase asymptotic',
                    'CP_Ma': 'Continuous Phase asymptotic'}
    if table_name in name_dict:
        return name_dict[table_name]
    else: return table_name
def output_graphs(kwargs):
    cur, conn = connect(kwargs['conf'],kwargs)
    # test_items = generate_testitems(kwargs)
    cur.execute(
        " SELECT name from sqlite_master WHERE type='table'")
    table_names = cur.fetchall()
    for wrapped_table_name in table_names:
        table_name = wrapped_table_name[0]
        if table_name in discards: discard = discards[table_name]
        else: discard = {}
        draw_single_line(
            cur, table_name, graph_name(table_name), discard)
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
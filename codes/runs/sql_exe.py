from omegaconf import OmegaConf
from sqlite3 import Cursor
import sqlite3
from root import ROOT_PATH
import json
import numpy as np

def table_exists(cur: Cursor, table_name):
    cur.execute(
        " SELECT name from sqlite_master WHERE type='table' AND name=?", (table_name,))
    res = cur.fetchall()
    if res == []:
        return False
    else:
        return True


def create_table(cur: Cursor,table_name, kwargs):
    if not table_exists(cur, table_name):
        cur.execute('''CREATE TABLE %s
                    (
                    x REAL[],
                    length REAL,
                    keyrate REAL,
                    CONSTRAINT constraint_name PRIMARY KEY (length)
                    );'''%table_name)


def search(cur: Cursor, table_name, len=None):
    if len != None:
        cur.execute('''SELECT * FROM '%s' WHERE length = ?''' %
                    table_name, (len,))
    else:
        cur.execute('SELECT * FROM "%s" ORDER BY length ASC'%table_name)
    res = cur.fetchall()
    return res

def x_to_sql_x(x):
    x= np.array(x)
    return '{'+json.dumps(x.tolist())[1:-1]+'}'
def write_records(cur: Cursor, x, len, keyrate, table_name):
    is_exist = search(cur, table_name, len)
    if is_exist != []:
        if is_exist[0][2] < keyrate:
            cur.execute('''DELETE FROM '%s'
                        where length = ? 
                        ''' % table_name, (len,))

            cur.execute('''INSERT INTO '%s' (x,length,keyrate)
            VALUES (?,?,?)'''%table_name, (x_to_sql_x(x), len, keyrate))
    else:
        cur.execute('''INSERT INTO '%s'(x,length,keyrate)
        VALUES (?,?,?)''' % table_name, (x_to_sql_x(x), len, keyrate))


def connect(conf,kwargs):
    base_name = ROOT_PATH+'/databases/'+conf+'.db'
    conn = sqlite3.connect(base_name)
    cur = conn.cursor()
    return cur,conn

def get_initials():
    base_name = ROOT_PATH+'/databases/initials.db'
    conn = sqlite3.connect(base_name)
    cur = conn.cursor()
    table_existed = table_exists(cur,'initials')
    if(table_existed):
        cur.execute('''SELECT * FROM initials''' )
        res = cur.fetchall()
    else:
        cur.execute(''' CREATE TABLE initials(
                        table_name TEXT,
                        x REAL[]
                        ); 
                    ''' )
        conn.commit()
        conn.close()
        res = []
    return res
def put_initials(table_name, x):
    base_name = ROOT_PATH+'/databases/initials.db'
    conn = sqlite3.connect(base_name)
    cur = conn.cursor()
    table_existed = table_exists(cur,'initials')
    if(not table_existed):
       cur.execute(''' CREATE TABLE initials(
                        table_name TEXT,
                        x REAL[]
                        ); 
                    ''' )
    cur.execute(
    " SELECT * from initials WHERE table_name=?", (table_name,))
    temp = cur.fetchall()
    if(temp!=[]):
        cur.execute(
        " UPDATE initials SET x = ? WHERE table_name=?", (x_to_sql_x(x),table_name,))
    elif x!=[]:
        cur.execute(
        " INSERT INTO initials VALUES (?,?)", (table_name,x_to_sql_x(x),))
    conn.commit()
    conn.close()


def get_table_name(cur):
    cur.execute('''SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;''')
    res = cur.fetchall()
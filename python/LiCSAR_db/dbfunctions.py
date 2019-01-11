#!/usr/bin/env python

import pymysql
from ConfigParser import SafeConfigParser
import global_config as gc


def Conn_db():
    # Returns the connection handle
    # Get DB connection info from config file
    parser = SafeConfigParser()
    parser.read(gc.configfile)
    sqlhost = parser.get('sqlinfo', 'host')
    sqldb = parser.get('sqlinfo', 'dbname')
    sqluser = parser.get('sqlinfo', 'dbuser')
    sqlpass = parser.get('sqlinfo', 'dbpass')

    try:
        # Connect to the database
        conn = pymysql.connect(host=sqlhost,
                               user=sqluser,
                               password=sqlpass,
                               db=sqldb)

        cur = conn.cursor()
        cur.execute('SELECT VERSION();')
        res = cur.fetchone()
        if not res:
            print "Error in database communication"
            return 'MYSQL ERROR'

    except pymysql.err.OperationalError as E:
        print 'MySQL ERROR:\n%s\n' % E
        conn = 'MYSQL ERROR'

    return conn

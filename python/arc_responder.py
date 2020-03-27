#!/usr/bin/env python

# Mar 2020 - Milan Lazecky

import sshtunnel
import pymysql
import pandas as pd

# first step to run through ssh tunnel to CEDA, in order to connect to LiCSInfo
#just do not forget about security J
sql_hostname = '130.246.129.18'
sql_username = 'lics'
sql_password = 'T34mLiCS'
sql_batch_database = 'licsinfo_batch'
sql_live_database = 'licsinfo_live'
sql_port = 3306

ssh_host = 'cems-login1.cems.rl.ac.uk'
ssh_username = 'earmla'
ssh_pkey = '~/.ssh/id_jasmin'
ssh_port = 22

tunnel = sshtunnel.SSHTunnelForwarder((ssh_host, ssh_port),
            ssh_username=ssh_username, ssh_pkey=ssh_pkey,
            remote_bind_address=(sql_hostname, sql_port))

#tunnel.start()  #this works, but in case of error it may be kept ON?
with tunnel:
    conn_batch = pymysql.connect(host='127.0.0.1',
                 user=sql_username,
                 passwd=sql_password,
                 database=sql_batch_database,
                 port=tunnel.local_bind_port)
    data = pd.read_sql_query("SHOW DATABASES;", conn_batch)

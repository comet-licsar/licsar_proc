#!/bin/bash
cd /gws/smf/j04/nceo_geohazards/software/iri2020/src/iri2020/data
mv apf107.dat apf107.dat.backup
wget https://chain-new.chain-project.net/echaim_downloads/apf107.dat  # daily update
#wget http://irimodel.org/indices/apf107.dat  # twice a year update
mv ig_rz.dat ig_rz.dat.backup
wget https://chain-new.chain-project.net/echaim_downloads/ig_rz.dat
#wget http://irimodel.org/indices/ig_rz.dat

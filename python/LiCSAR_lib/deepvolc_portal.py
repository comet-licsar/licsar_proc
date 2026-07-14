#!/usr/bin/env python
'''
- first get the nc file:
LiCSBAS_out2nc.py -i TS_GEOCml2clipmask/cum.h5 --cf -o $LBDIR/deepvolc_portal/test_083D.nc
- then copy to the server:
scp -A  $LBDIR/deepvolc_portal/test_083D.nc rocky@130.246.213.180:/mnt/deepvolc_1/test1/.
- then inside the server, ingest it to ncWMS:
dataid='test_083D'
adus=licsadmin
adpas=`cat /home/rocky/.nc.ps`
curl --digest -u $adus:$adpas \
  -d "id="$dataid"&title="$dataid"&location=/mnt/deepvolc_1/test1/"$dataid".nc" \
  -X POST http://localhost:8080/ncWMS2/admin/addDataset

# make it downloadable:
curl --digest -u $adus:$adpas \
  -X POST \
  -d "id=$dataid" \
  -d "downloadable=true" \
  -d "queryable=true" \
  http://localhost:8080/ncWMS2/admin/updateDataset


it works? check:
http://130.246.213.180/ncWMS2/

also it now should be downloadable from:
http://130.246.213.180/deepvolc/test_083D.nc

now, nginx expects XCUBE listening on 8443 (finally!!! see /etc/nginx/..., so I can try start it using:
 micromamba activate xcube
 xcube serve --verbose --port 8443 /mnt/deepvolc_1/test1/test_tmp.nc

now xcube should be reachable from
https://130.246.213.180.sslip.io/viewer

'''
#!/bin/bash
set -x


branch=$(git branch -a |grep '*' | awk -F ' ' '{print $2}')


#RCP="benshila@legos.obs-mip.fr,patrick.marchesiello@ird.fr,gildas.cambon@ird.fr"
RCP="benshila@legos.obs-mip.fr,patrick.marchesiello@ird.fr,gildas.cambon@ird.fr"

cat <<EOF > msg.txt
Tests results for $(date +%Y-%m-%d) (right) compared to $(date -d "-1 days" +%Y-%m-%d) (left) for branch ${branch}.

This an automated message. DO NOT ANSWER. 


EOF


#mutt -e 'my_hdr From: CROCO tests <croco_tests@do_not_reply>' -s "CROCO tests $(date +%Y-%m-%d)" -i msg.txt -a ./PDF/diffplot_$(date +%Y-%m-%d).pdf  -- $RCP < /dev/null
mutt -e 'set realname ="CROCO test" ' -e 'set from=DoNotReply' -s "CROCO tests for $branch $(date +%Y-%m-%d)" -i msg.txt -a ./PDF/diffplot_$(date +%Y-%m-%d).pdf  -- $RCP < /dev/null

\rm msg.txt

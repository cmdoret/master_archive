# This is a simple script to look
# Cyril Matthey-Doret
# 20.10.2017

orig_ref=
order_ref=


# Convert ordered assembly coordinates to original assembly and look up
# annotations around the coordinate
for hit in ;
do
  ori_hit=$(python2 chr2contig.py $hit $order_ref $orig_ref --region_size 10000)
  python2 pos2annot.py $ori_hit
done

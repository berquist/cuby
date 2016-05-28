c=0 
grep " name:.*1.00" $1 | sed -e 's/^ *name: //' -e 's/ 1.00.*//' | while read i
do 
	f=$(($c*8))
	l=$(($f+7))
	fn=`echo $i | sed -e 's/ \.\.\. /-/' -e 's/ /_/'`
	sed -e "s/%name%/$i/" -e "s/%first%/$f/" -e "s/%last%/$l/" -e "s/%fn%/$fn/" source_tables/curve_template8
       	c=$(($c+1))
done

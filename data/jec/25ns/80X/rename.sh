for f in *.txt; do 
    cp "$f" "`echo $f | sed s/Spring16_25nsV6_//`"; 
done
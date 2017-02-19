for f in *.txt; do 
    cp "$f" "`echo $f | sed s/Summer16_23Sep2016V4_//`"; 
done

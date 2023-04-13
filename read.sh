#!/bin/bash
input=$1
while IFS= read -r line
 do
 #echo "$line"
root4star -l -b <<EOF
.x HelloWorld.C(1,("$line"))

.q
EOF
done <$input
	   

#!/bin/bash

for((i=1;i<=10;i++))
do
	dos2unix seq${i}.fa
	sed -i "s/^M//" seq${i}.fa
done
echo done!

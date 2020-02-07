#!/bin/bash -x

#for d in vs* epsc* pred* ssteps* minsteps* spred* ; do 
for d in ssteps[23456789]/ ; do 
	cp -p runFinLearn.R utils.R worlds.R FinanceDemo.R finlearn.R findWayFinLearn2.R $d/; 
done

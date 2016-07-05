#!/bin/bash

# Ryan Y
# Script takes as an argument a folder or path ($1) and synchronizes them on this
# and a remote computer

rsync -av ryoung@129.64.45.231:~/Data/Miller/$1 $MData
rsync -av $MData/$1 ryoung@129.64.45.231:~/Data/Miller/

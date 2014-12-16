#!/bin/bash

p1111=`ls $1 | grep '1-1--1-1'` 
p2211=`ls $1 | grep '2-2--1-1'` 
p1122=`ls $1 | grep '1-1--2-2'` 
p2222=`ls $1 | grep '2-2--2-2'` 
p1212=`ls $1 | grep '1-2--1-2'`

p1112=`ls $1 | grep '1-1--1-2'`
p2212=`ls $1 | grep '2-2--1-2'`
p1211=`ls $1 | grep '1-2--1-1'`
p1222=`ls $1 | grep '1-2--2-2'`


# secular contributions
xmgrace $p1111 $p2211 $p1122 $p2222 $p1212 &

# non-secular contributions
xmgrace $p1112 $p2212 &

#!/bin/bash

wget -U Mozilla http://automation.whatismyip.com/n09230945.asp
echo -n "mini: "   > $HOME/dropbox/ip.txt
cat n09230945.asp >> $HOME/dropbox/ip.txt 
rm -f n09230945.asp 
scp $HOME/dropbox/ip.txt solidpress@solidstatepress.com:html/ip.html
